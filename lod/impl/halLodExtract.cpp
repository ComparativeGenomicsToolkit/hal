/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <deque>
#include <limits>
#include <algorithm>
#include "halLodExtract.h"
extern "C" {
#include "sonLibTree.h"
}
 

using namespace std;
using namespace hal;

LodExtract::LodExtract()
{

}

LodExtract::~LodExtract()
{
  
}

void LodExtract::createInterpolatedAlignment(AlignmentConstPtr inAlignment,
                                             AlignmentPtr outAlignment,
                                             double scale,
                                             const string& tree,
                                             const string& rootName,
                                             bool keepSequences,
                                             bool allSequences,
                                             double probeFrac,
                                             double minSeqFrac)
{
  _inAlignment = inAlignment;
  _outAlignment = outAlignment;
  _keepSequences = keepSequences;
  _allSequences = allSequences;
  _probeFrac = probeFrac;
  _minSeqFrac = minSeqFrac;
  
  string newTree = tree.empty() ? inAlignment->getNewickTree() : tree;
  createTree(newTree, rootName);
  cout << "tree = " << _outAlignment->getNewickTree() << endl;
  
  deque<string> bfQueue;
  bfQueue.push_front(_outAlignment->getRootName());
  while (!bfQueue.empty())
  {
    string genomeName = bfQueue.back();
    bfQueue.pop_back();
    vector<string> childNames = _outAlignment->getChildNames(genomeName);
    if (!childNames.empty())
    {
      convertInternalNode(genomeName, scale);
      for (size_t childIdx = 0; childIdx < childNames.size(); childIdx++)
      {
        bfQueue.push_back(childNames[childIdx]);
      } 
    }
  }
}

void LodExtract::createTree(const string& tree, const string& rootName)
{
  if (_outAlignment->getNumGenomes() != 0)
  {
    throw hal_exception("Output alignment not empty");
  }
  stTree* treeRoot = stTree_parseNewickString(const_cast<char*>(tree.c_str()));
  stTree* root = treeRoot;
  assert(treeRoot != NULL);
  deque<stTree*> bfQueue;
  
  // find rootName in tree if specified, otherwise we use the whole tree
  if (rootName.empty() == false)
  {
    bfQueue.push_front(treeRoot);
    while (!bfQueue.empty())
    {
      stTree* node = bfQueue.back();
      const char* label = stTree_getLabel(node);
      if (rootName == string(label))
      {
        root = node;
        break;
      }
      int32_t numChildren = stTree_getChildNumber(node);
      for (int32_t childIdx = 0; childIdx < numChildren; ++childIdx)
      {
        bfQueue.push_front(stTree_getChild(node, childIdx));
      }

      bfQueue.pop_back();
    }
    bfQueue.clear();
  }

  bfQueue.push_front(root);
  while (!bfQueue.empty())
  {
    stTree* node = bfQueue.back();
    const char* label = stTree_getLabel(node);
    if (label == NULL)
    {
      throw hal_exception("Error parsing tree: unlabeled node");
    }
    const Genome* test = _inAlignment->openGenome(label);
    if (test == NULL)
    {
      throw hal_exception(string("Genome in tree: ") + string(label) +
                          "doesn't exist in source alignment");
    }
    else
    {
      _inAlignment->closeGenome(test);
    }    
    if (node == root)
    {
      _outAlignment->addRootGenome(label);
    }
    else
    {
      stTree* parent = stTree_getParent(node);
      assert(parent != NULL);
      const char* pLabel = stTree_getLabel(parent);
      assert(pLabel != NULL);
      double branchLength = stTree_getBranchLength(node);
      // clamp undefined branch lengths to 1. for now
      if (branchLength > 1e10)
      {
        branchLength = 1.;
      }
      _outAlignment->addLeafGenome(label, pLabel, branchLength);
    }

    int32_t numChildren = stTree_getChildNumber(node);
    for (int32_t childIdx = 0; childIdx < numChildren; ++childIdx)
    {
      bfQueue.push_front(stTree_getChild(node, childIdx));
    }

    bfQueue.pop_back();
  }
  stTree_destruct(root);
}

void LodExtract::convertInternalNode(const string& genomeName, 
                                     double scale)
{
  const Genome* parent = _inAlignment->openGenome(genomeName);
  assert(parent != NULL);
  vector<string> childNames = _outAlignment->getChildNames(genomeName);
  vector<const Genome*> children;
  for (hal_size_t i = 0; i < childNames.size(); ++i)
  {
    children.push_back(_inAlignment->openGenome(childNames[i]));
  }
  const Genome* grandParent = NULL; // TEMP HACK  parent->getParent();
  hal_size_t minAvgBlockSize = getMinAvgBlockSize(parent, children, grandParent);
  hal_size_t step = (hal_size_t)(scale * minAvgBlockSize);
  _graph.build(_inAlignment, parent, children, grandParent, step, _allSequences, 
               _probeFrac, _minSeqFrac);

  map<const Sequence*, hal_size_t> segmentCounts;
  countSegmentsInGraph(segmentCounts);

  writeDimensions(segmentCounts, parent->getName(), childNames);
  if (_keepSequences == true)
  {
    writeSequences(parent, children);
  }
  writeSegments(parent, children);
  writeHomologies(parent, children);
  writeParseInfo(_outAlignment->openGenome(parent->getName()));

  // if we're gonna print anything out, do it before this:
  // (not necesssary but by closing genomes we erase their hdf5 caches
  // which can make a difference on huge trees
  _graph.erase();
  _outAlignment->closeGenome(_outAlignment->openGenome(parent->getName()));
  _inAlignment->closeGenome(parent);
  for (hal_size_t i = 0; i < children.size(); ++i)
  {
    _outAlignment->closeGenome(
      _outAlignment->openGenome(children[i]->getName()));
    _inAlignment->closeGenome(children[i]);    
  }
  if (grandParent != NULL)
  {
    _outAlignment->closeGenome(
      _outAlignment->openGenome(grandParent->getName()));
    _inAlignment->closeGenome(grandParent);
  }
}

void LodExtract::countSegmentsInGraph(
  map<const Sequence*, hal_size_t>& segmentCounts)
{
  const LodBlock* block;
  const LodSegment* segment;
  pair<map<const Sequence*, hal_size_t>::iterator, bool> res;

  for (hal_size_t blockIdx = 0; blockIdx < _graph.getNumBlocks(); ++blockIdx)
  {
    block = _graph.getBlock(blockIdx);
    for (hal_size_t segIdx = 0; segIdx < block->getNumSegments(); ++segIdx)
    {
      segment = block->getSegment(segIdx);
      res = segmentCounts.insert(pair<const Sequence*, hal_size_t>(
                                   segment->getSequence(), 0));
      hal_size_t& count = res.first->second;
      ++count;
    }
  }
  
  // add unsampled non-zero sequences to dimensions, by looking for
  // sequences who have telomeres but no segments. 
  const LodBlock* telomeres = _graph.getTelomeres();
  for (hal_size_t telIdx = 0; telIdx < telomeres->getNumSegments(); ++telIdx)
  {
    segment = telomeres->getSegment(telIdx);
    if (segment->getSequence()->getSequenceLength() > 0)
    {
      res = segmentCounts.insert(pair<const Sequence*, hal_size_t>(
                                   segment->getSequence(), 0));
      hal_size_t& count = res.first->second;
      if (res.second == true)
      {
        ++count;
        assert(count == 1);
      }
    }
  }
}

void LodExtract::writeDimensions(
  const map<const Sequence*, hal_size_t>& segmentCounts, 
  const string& parentName,
  const vector<string>& childNames)
{
  // initialize a dimensions list for each (input) genome
  map<const Genome*, vector<Sequence::Info> > dimMap;
  map<const Genome*, vector<Sequence::Info> >::iterator dimMapIt;
  vector<string> newGenomeNames = childNames;
  newGenomeNames.push_back(parentName);
 
  for (size_t i = 0; i < newGenomeNames.size(); ++i)
  {
    const Genome* inGenome = _inAlignment->openGenome(newGenomeNames[i]);
    pair<const Genome*, vector<Sequence::Info> > newEntry;
    newEntry.first = inGenome;
    
    // it's important we keep the sequences in the output genome
    // in the same order as the sequences in the input genome since
    // we always use global coordinates!
    SequenceIteratorConstPtr seqIt = inGenome->getSequenceIterator();
    SequenceIteratorConstPtr seqEnd = inGenome->getSequenceEndIterator();
    for (; seqIt != seqEnd; seqIt->toNext())
    {
      const Sequence* inSequence = seqIt->getSequence();
      map<const Sequence*, hal_size_t>::const_iterator segMapIt;
      segMapIt = segmentCounts.find(inSequence);
      // we skip empty sequences for now with below check
      if (segMapIt != segmentCounts.end())
      {
        vector<Sequence::Info>& segDims = newEntry.second;
        hal_size_t nTop = 
           inGenome->getName() == parentName ? 0 : segMapIt->second;
        hal_size_t nBot = 
           inGenome->getName() != parentName ? 0 : segMapIt->second;
        segDims.push_back(Sequence::Info(inSequence->getName(),
                                         inSequence->getSequenceLength(),
                                         nTop,
                                         nBot));
      }
    }

    // note potential bug here for genome with no data
    dimMap.insert(newEntry);
  }
  
  // now that we have the dimensions for each genome, update them in
  // the output alignment
  for (dimMapIt = dimMap.begin(); dimMapIt != dimMap.end(); ++dimMapIt)
  {
    Genome* newGenome = _outAlignment->openGenome(dimMapIt->first->getName());
    assert(newGenome != NULL);
    vector<Sequence::Info>& segDims = dimMapIt->second;
    // ROOT 
    if (newGenome->getName() == _outAlignment->getRootName())
    {
      assert(newGenome->getName() == parentName);
      newGenome->setDimensions(segDims, _keepSequences);
    }
    // LEAF
    else if (newGenome->getName() != parentName)
    {
      newGenome->setDimensions(segDims, _keepSequences);
    }
    // INTERNAL NODE
    else
    {
      vector<Sequence::UpdateInfo> updateInfo;
      for (size_t i = 0; i < segDims.size(); ++i)
      {
        updateInfo.push_back(
          Sequence::UpdateInfo(segDims[i]._name,
                               segDims[i]._numBottomSegments));
      }
      newGenome->updateBottomDimensions(updateInfo);
    }
  }
}

void LodExtract::writeSequences(const Genome* inParent,
                                const vector<const Genome*>& inChildren)
{
  vector<const Genome*> inGenomes = inChildren;
  inGenomes.push_back(inParent);
  const Genome* outParent = _outAlignment->openGenome(inParent->getName());
  (void)outParent;
  assert(outParent != NULL && outParent->getNumBottomSegments() > 0);
  string buffer;

  for (hal_size_t i = 0; i < inGenomes.size(); ++i)
  {
    const Genome* inGenome = inGenomes[i];
    Genome* outGenome = _outAlignment->openGenome(inGenome->getName());
    if (inGenome == inParent || outGenome->getNumChildren() == 0)
    {   
      SequenceIteratorConstPtr inSeqIt = inGenome->getSequenceIterator();
      SequenceIteratorConstPtr end = inGenome->getSequenceEndIterator();
      for (; inSeqIt != end; inSeqIt->toNext())
      {
        const Sequence* inSequence = inSeqIt->getSequence();
        if (inSequence->getSequenceLength() > 0)
        {
          Sequence* outSequence = outGenome->getSequence(inSequence->getName());
          assert(outSequence != NULL);
          inSequence->getString(buffer);
          outSequence->setString(buffer);
        }
      }
    }
  }
}

void LodExtract::writeSegments(const Genome* inParent,
                               const vector<const Genome*>& inChildren)
{
  vector<const Genome*> inGenomes = inChildren;
  inGenomes.push_back(inParent);
  const Genome* outParent = _outAlignment->openGenome(inParent->getName());
  assert(outParent != NULL && outParent->getNumBottomSegments() > 0);
  BottomSegmentIteratorPtr bottom;
  TopSegmentIteratorPtr top;
  SegmentIteratorPtr outSegment;

  // FOR EVERY GENOME
  for (hal_size_t i = 0; i < inGenomes.size(); ++i)
  {
    const Genome* inGenome = inGenomes[i];
    Genome* outGenome = _outAlignment->openGenome(inGenome->getName());

    SequenceIteratorPtr outSeqIt = outGenome->getSequenceIterator();
    SequenceIteratorConstPtr outSeqEnd = outGenome->getSequenceEndIterator();
    
    // FOR EVERY SEQUENCE IN GENOME
    for (; outSeqIt != outSeqEnd; outSeqIt->toNext())
    {
      const Sequence* outSequence = outSeqIt->getSequence();
      const Sequence* inSequence = 
         inGenome->getSequence(outSequence->getName());
      if (outGenome != outParent && outSequence->getNumTopSegments() > 0)
      {
        top = outSequence->getTopSegmentIterator();
        outSegment = top;
      }
      else if (outSequence->getNumBottomSegments() > 0)
      {
        bottom = outSequence->getBottomSegmentIterator();
        outSegment = bottom;
      }
      const LodGraph::SegmentSet* segSet = _graph.getSegmentSet(inSequence);
      assert(segSet != NULL);
      LodGraph::SegmentSet::const_iterator segIt = segSet->begin();
      if (segSet->size() > 2)
      {
        //skip left telomere
        ++segIt;
        // use to skip right telomere:
        LodGraph::SegmentSet::const_iterator segLast = segSet->end();
        --segLast;
      
        // FOR EVERY SEGMENT IN SEQUENCE
        for (; segIt != segLast; ++segIt)
        {
          // write the HAL array index back to the segment to make
          // future passes quicker. 
          (*segIt)->setArrayIndex(outSegment->getArrayIndex());
          outSegment->setCoordinates((*segIt)->getLeftPos(), 
                                     (*segIt)->getLength());
          assert(outSegment->getSequence()->getName() == inSequence->getName());
          outSegment->toRight();
        }
      }
      else if (outSequence->getSequenceLength() > 0)
      {
        assert(segSet->size() == 2);
        writeUnsampledSequence(outSequence, outSegment);
      }
    }
  } 
}

void LodExtract::writeUnsampledSequence(const Sequence* outSequence,
                                        SegmentIteratorPtr outSegment)
{
  outSegment->setCoordinates(outSequence->getStartPosition(),
                             outSequence->getSequenceLength());
  if (outSegment->isTop())
  {
    assert(outSequence->getNumTopSegments() == 1);
    TopSegmentIteratorPtr top = outSegment.downCast<TopSegmentIteratorPtr>();
    top->setParentIndex(NULL_INDEX);
    top->setParentReversed(false);
    top->setNextParalogyIndex(NULL_INDEX);
    top->setBottomParseIndex(NULL_INDEX);
  }
  else
  {
    assert(outSequence->getNumBottomSegments() == 1);
    BottomSegmentIteratorPtr bottom = 
       outSegment.downCast<BottomSegmentIteratorPtr>();
    hal_size_t numChildren = bottom->getNumChildren();
    for (hal_size_t childNum = 0; childNum < numChildren; ++childNum)
    {
      bottom->setChildIndex(childNum, NULL_INDEX);
      bottom->setChildReversed(childNum, false);
    }
    bottom->setTopParseIndex(NULL_INDEX);
  }
}

void LodExtract::writeHomologies(const Genome* inParent,
                                 const vector<const Genome*>& inChildren)
{
  vector<const Genome*> inGenomes = inChildren;
  inGenomes.push_back(inParent);
  Genome* outParent = _outAlignment->openGenome(inParent->getName());
  assert(outParent != NULL && outParent->getNumBottomSegments() > 0);
  assert(inChildren.size() > 0);
  Genome* outChild = _outAlignment->openGenome(inChildren[0]->getName());
  BottomSegmentIteratorPtr bottom = outParent->getBottomSegmentIterator();
  TopSegmentIteratorPtr top = outChild->getTopSegmentIterator();

  // FOR EVERY BLOCK
  for (hal_size_t blockIdx = 0; blockIdx < _graph.getNumBlocks(); ++blockIdx)
  {
    SegmentMap segMap;
    const LodBlock* block = _graph.getBlock(blockIdx);

    for (hal_size_t segIdx = 0; segIdx < block->getNumSegments(); ++segIdx)
    {
      const LodSegment* segment = block->getSegment(segIdx);
      const Genome* genome = segment->getSequence()->getGenome();

      // ADD TO MAP
      pair<SegmentMap::iterator, bool> res = segMap.insert(
        pair<const Genome*, SegmentSet*>(genome, NULL));
      if (res.second == true)
      {
        assert(res.first->second == NULL);
        res.first->second = new SegmentSet();
      }
      res.first->second->insert(segment);    
    }      
    updateBlockEdges(inParent, segMap, block, bottom, top);
    
    // free the temporary sets! 
    for (SegmentMap::iterator mapIt = segMap.begin(); mapIt != segMap.end();
         ++mapIt)
    {
      delete mapIt->second;
    }
  }
}

void LodExtract::updateBlockEdges(const Genome* inParentGenome,
                                  SegmentMap& segMap,
                                  const LodBlock* block,
                                  BottomSegmentIteratorPtr bottom,
                                  TopSegmentIteratorPtr top)
{
  Genome* outParentGenome = bottom->getGenome();
  const LodSegment* rootSeg = NULL;
  SegmentSet* segSet;
  SegmentSet::iterator setIt;

  // Zap all segments in parent genome
  SegmentMap::iterator mapIt = segMap.find(inParentGenome);
  if (mapIt != segMap.end())
  {
    segSet = mapIt->second;
    assert(segSet != NULL);
    setIt = segSet->begin();
    for (; setIt != segSet->end(); ++setIt)
    {
      bottom->setArrayIndex(outParentGenome, (*setIt)->getArrayIndex());
      for (hal_size_t i = 0; i < bottom->getNumChildren(); ++i)
      {
        bottom->setChildIndex(i, NULL_INDEX);
        bottom->setTopParseIndex(NULL_INDEX);
      }
    }

    // Choose first segment as parent to all segments in the child genome
    setIt = segSet->begin();
    rootSeg = *(setIt);
    bottom->setArrayIndex(outParentGenome, (*setIt)->getArrayIndex());
  }

  // Do the child genomes
  const Genome* inGrandParentGenome = inParentGenome->getParent();
  SegmentSet::iterator nextIt;
  for (mapIt = segMap.begin(); mapIt != segMap.end(); ++mapIt)
  {
    if (mapIt->first != inParentGenome and mapIt->first != inGrandParentGenome)
    {
      Genome* outChildGenome = 
         _outAlignment->openGenome(mapIt->first->getName());
      hal_index_t childIndex = outParentGenome->getChildIndex(outChildGenome);
      assert(childIndex >= 0);
      segSet = mapIt->second;
      assert(segSet != NULL);
      for (setIt = segSet->begin(); setIt != segSet->end(); ++setIt)
      {
        top->setArrayIndex(outChildGenome, (*setIt)->getArrayIndex());
        top->setBottomParseIndex(NULL_INDEX);

        // Connect to parent
        if (rootSeg != NULL)
        {
          top->setParentIndex(bottom->getArrayIndex());
          bool reversed = (*setIt)->getFlipped() != rootSeg->getFlipped();
          top->setParentReversed(reversed);
          if (setIt == segSet->begin())
          {
            bottom->setChildIndex(childIndex, top->getArrayIndex());         
            bottom->setChildReversed(childIndex, reversed);      
          }
        }
        else
        {
          top->setParentIndex(NULL_INDEX);
        }

        // Connect to next paralogy
        SegmentSet::iterator setNext = setIt;
        ++setNext;
        if (setNext == segSet->end())
        {
          setNext = segSet->begin();
        }
        if (setNext == setIt)
        {
          top->setNextParalogyIndex(NULL_INDEX);
        }
        else
        {
          top->setNextParalogyIndex((*setNext)->getArrayIndex());
        }
      }
    }
  }
}

void LodExtract::writeParseInfo(Genome* genome)
{
  if (genome->getParent() == NULL || genome->getNumChildren() == 0)
  {
    return;
  }

 // copied from CactusHalConverter::updateRootParseInfo() in
  // cactus2hal/src/cactusHalConverter.cpp 
  BottomSegmentIteratorPtr bottomIterator = 
     genome->getBottomSegmentIterator();
  TopSegmentIteratorPtr topIterator = genome->getTopSegmentIterator();
  BottomSegmentIteratorConstPtr bend = genome->getBottomSegmentEndIterator();
  TopSegmentIteratorConstPtr tend = genome->getTopSegmentEndIterator();

  while (bottomIterator != bend && topIterator != tend)
  {
    bool bright = false;
    bool tright = false;
    BottomSegment* bseg = bottomIterator->getBottomSegment();
    TopSegment* tseg = topIterator->getTopSegment();
    hal_index_t bstart = bseg->getStartPosition();
    hal_index_t bend = bstart + (hal_index_t)bseg->getLength();
    hal_index_t tstart = tseg->getStartPosition();
    hal_index_t tend = tstart + (hal_index_t)tseg->getLength();
    
    if (bstart >= tstart && bstart < tend)
    {
      bseg->setTopParseIndex(tseg->getArrayIndex());
    }
    if (bend <= tend || bstart == bend)
    {
      bright = true;
    }
        
    if (tstart >= bstart && tstart < bend)
    {
      tseg->setBottomParseIndex(bseg->getArrayIndex());
    }
    if (tend <= bend || tstart == tend)
    {
      tright = true;
    }

    assert(bright || tright);
    if (bright == true)
    {
      bottomIterator->toRight();
    }
    if (tright == true)
    {
      topIterator->toRight();
    }
  }
}


hal_size_t LodExtract::getMinAvgBlockSize(
  const Genome* inParent,
  const vector<const Genome*>& inChildren,
  const Genome* inGrandParent) const
{
  hal_size_t minAvgBlockSize = numeric_limits<hal_size_t>::max();
  if (inParent->getSequenceLength() > 0)
  {
    assert(inParent->getNumBottomSegments() > 0);
    assert(inParent->getNumBottomSegments() < inParent->getSequenceLength());
    minAvgBlockSize = inParent->getSequenceLength() / 
       inParent->getNumBottomSegments();
  }
  for (size_t i = 0; i < inChildren.size(); ++i)
  {
    if (inChildren[i]->getSequenceLength() > 0)
    {
      assert(inChildren[i]->getNumTopSegments() > 0);
      assert(inChildren[i]->getNumTopSegments() <
             inParent->getSequenceLength());
      minAvgBlockSize = std::min(minAvgBlockSize, 
                                 inChildren[i]->getSequenceLength() / 
                                 inChildren[i]->getNumTopSegments());
    }
  }
  if (inGrandParent != NULL && inGrandParent->getSequenceLength() > 0)
  {
    assert (inGrandParent->getNumBottomSegments() > 0);
    assert (inGrandParent->getNumBottomSegments() < 
            inGrandParent->getSequenceLength());
    minAvgBlockSize = std::min(minAvgBlockSize,
                               inGrandParent->getSequenceLength() /
                               inGrandParent->getNumBottomSegments());
  }
  return minAvgBlockSize;
}
