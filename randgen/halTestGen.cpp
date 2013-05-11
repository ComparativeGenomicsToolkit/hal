/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include "hal.h"
#include "halTopSegmentTest.h"
#include "halBottomSegmentTest.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "output hal file");
  return optionsParser;
}

static void setDNA(Genome* genome, char base = 'g');
static void writeParseInfo(Genome* genome);
static void createRoot(AlignmentPtr alignment, const string& name,
                       hal_size_t numBot, hal_size_t botWidth);
static void addLeaf(AlignmentPtr alignment, const string& name,
                    hal_size_t numTop);
static void invertLeaf(AlignmentPtr alignment, 
                       hal_index_t firstTop, hal_index_t lastTop);
static void dupLeaf(AlignmentPtr alignment,
                    hal_size_t dupLen, hal_index_t sourceIdx,
                    const vector<hal_index_t>& tgtIdx,
                    const vector<bool>& tgtRev);
static void transLeaf(AlignmentPtr alignment,
                      hal_size_t len, hal_index_t srcIdx,
                      hal_index_t tgtIdx);
static void indelLeaf(AlignmentPtr alignment,
                      hal_index_t firstTop, hal_index_t lastTop);

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string halPath;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  try
  {
    const hal_size_t genomeLength = 5000;
    hal_size_t segSize;
    
    AlignmentPtr alignment = hdf5AlignmentInstance();
    alignment->setOptionsFromParser(optionsParser);
    alignment->createNew(halPath);
    
    createRoot(alignment, "root", 1, genomeLength);
    
    segSize = 10;
    addLeaf(alignment, "1inv_1000_1500", genomeLength / segSize);
    invertLeaf(alignment, 1000 / segSize, (1500 - segSize) / segSize);

    segSize = 20;
    addLeaf(alignment, "2inv_1200_1400", genomeLength / segSize);
    invertLeaf(alignment, 1200 / segSize, (1400 - segSize) / segSize);

    segSize = 100;
    addLeaf(alignment, "3dup_200_500_600_800r", genomeLength / segSize);
    vector<hal_index_t> tgtIdx(3);
    vector<bool> tgtRev(3, false);
    tgtIdx[0] = 500 / segSize;
    tgtIdx[1] = 600 / segSize;
    tgtIdx[2] = 800 /segSize;
    tgtRev[2] = true;
    dupLeaf(alignment, 1, 200 / segSize, tgtIdx, tgtRev);

    segSize = 5;
    addLeaf(alignment, "4swap_20_920to3000_3900", genomeLength / segSize);
    transLeaf(alignment, 900 / segSize, 20 / segSize, 3000 / segSize);

    segSize = 1000;
    addLeaf(alignment, "5gap_0_1000", genomeLength / segSize);
    indelLeaf(alignment, 0 / segSize, (1000 - segSize) / segSize);

    segSize = 500;
    addLeaf(alignment, "6inv_2500_5000", genomeLength / segSize);
    invertLeaf(alignment, 2500 / segSize, (5000 - segSize) / segSize);

    validateAlignment(alignment);
    alignment->close();
  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  
  return 0;
}

void setDNA(Genome* genome, char base)
{
  string dna(genome->getSequenceLength(), base);
  genome->setString(dna);
}

void createRoot(AlignmentPtr alignment, const string& name,
                hal_size_t numBot, hal_size_t botWidth)
{
  assert(alignment->getNumGenomes() == 0);
  Genome* root = alignment->addRootGenome(name);
  Sequence::Info dims(name, botWidth * numBot, 0, numBot);
  root->setDimensions(vector<Sequence::Info>(1, dims));
  setDNA(root);
}

void addLeaf(AlignmentPtr alignment, const string& name,
             hal_size_t numTop)
{
  Genome* lastLeaf = alignment->openGenome(alignment->getRootName());
  while (lastLeaf->getNumChildren() > 0) 
  {
    assert(lastLeaf->getSequenceLength() % numTop == 0);
    lastLeaf = lastLeaf->getChild(0);
  }
  Genome* leaf = alignment->addLeafGenome(name, lastLeaf->getName(), 1);
  
  Sequence::UpdateInfo udims(lastLeaf->getName(), numTop);
  lastLeaf->updateBottomDimensions(vector<Sequence::UpdateInfo>(1, udims));

  Sequence::Info dims(name, lastLeaf->getSequenceLength(), numTop, 0);
  leaf->setDimensions(vector<Sequence::Info>(1, dims));
  setDNA(leaf);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  hal_index_t segLen = leaf->getSequenceLength() / numTop;
  bi = lastLeaf->getBottomSegmentIterator();
  ti = leaf->getTopSegmentIterator();
  for (; bi != lastLeaf->getBottomSegmentEndIterator(); 
       bi->toRight(), ti->toRight())
  {
    bs.set(bi->getArrayIndex() * segLen, segLen);
    bs._children.clear();
    bs._children.push_back(pair<hal_size_t, bool>(bi->getArrayIndex(), false));
    bs.applyTo(bi);

    ts.set(ti->getArrayIndex() * segLen, segLen, ti->getArrayIndex());
    ts.applyTo(ti);
  }
  writeParseInfo(lastLeaf);
}

void invertLeaf(AlignmentPtr alignment, 
                hal_index_t firstTop, hal_index_t lastTop)
{
  Genome* leaf = alignment->openGenome(alignment->getRootName());
  while (leaf->getNumChildren() > 0) 
  {
    leaf = leaf->getChild(0);
  }
  Genome* parent = leaf->getParent();  

  assert(firstTop <= lastTop);
  assert(lastTop < leaf->getNumTopSegments());
  assert(lastTop < parent->getNumBottomSegments());

  hal_index_t span = 1 + lastTop - firstTop;
  for (hal_index_t i = 0; i < span; ++i)
  {
    TopSegmentIteratorPtr top = leaf->getTopSegmentIterator(firstTop + i);
    BottomSegmentIteratorPtr bot =
       parent->getBottomSegmentIterator(lastTop - i);
    top->setParentIndex(bot->getArrayIndex());
    top->setParentReversed(true);
    bot->setChildIndex(0, top->getArrayIndex());
    bot->setChildReversed(0, true);
  }
}

void dupLeaf(AlignmentPtr alignment,
             hal_size_t dupLen, hal_index_t sourceIdx,
             const vector<hal_index_t>& tgtIdx,
             const vector<bool>& tgtRev)
{
  Genome* leaf = alignment->openGenome(alignment->getRootName());
  while (leaf->getNumChildren() > 0) 
  {
    leaf = leaf->getChild(0);
  }
  Genome* parent = leaf->getParent();  

  assert(sourceIdx <= leaf->getNumTopSegments());
  assert(tgtIdx.size() == tgtRev.size());

  for (size_t i = 0; i < tgtIdx.size(); ++i)
  {
    assert(tgtIdx[i] != sourceIdx);
    TopSegmentIteratorPtr top = leaf->getTopSegmentIterator(tgtIdx[i]);
    hal_index_t prevIdx = i == 0 ? sourceIdx : tgtIdx[i-1];
    TopSegmentIteratorPtr prev = leaf->getTopSegmentIterator(prevIdx);
    BottomSegmentIteratorPtr bot =
       parent->getBottomSegmentIterator(top->getParentIndex());
    
    top->setParentIndex(prev->getParentIndex());
    top->setParentReversed(tgtRev[i]);
    bot->setChildIndex(0, NULL_INDEX);

    prev->setNextParalogyIndex(top->getArrayIndex());
    top->setNextParalogyIndex(NULL_INDEX);
    if (i == tgtIdx.size() - 1)
    {
      top->setNextParalogyIndex(sourceIdx);
    }
  }
}

void transLeaf(AlignmentPtr alignment,
               hal_size_t len, hal_index_t srcIdx,
               hal_index_t tgtIdx)
{
  Genome* leaf = alignment->openGenome(alignment->getRootName());
  while (leaf->getNumChildren() > 0) 
  {
    leaf = leaf->getChild(0);
  }
  Genome* parent = leaf->getParent();  

  assert(tgtIdx < leaf->getNumTopSegments());
  assert(tgtIdx > len);
  assert(srcIdx < tgtIdx - len);
  assert(tgtIdx + len < leaf->getNumTopSegments());

  for (size_t i = 0; i < len; ++i)
  {
    TopSegmentIteratorPtr top = leaf->getTopSegmentIterator(srcIdx + i);
    TopSegmentIteratorPtr ttop = leaf->getTopSegmentIterator(tgtIdx + i);
    BottomSegmentIteratorPtr bot = parent->getBottomSegmentIterator(srcIdx + i);
    BottomSegmentIteratorPtr tbot =
       parent->getBottomSegmentIterator(tgtIdx + i);
    
    top->setParentIndex(tbot->getArrayIndex());
    tbot->setChildIndex(0, top->getArrayIndex());
    ttop->setParentIndex(bot->getArrayIndex());
    bot->setChildIndex(0, ttop->getArrayIndex());
  }
}

void indelLeaf(AlignmentPtr alignment,
               hal_index_t firstTop, hal_index_t lastTop)
{
  Genome* leaf = alignment->openGenome(alignment->getRootName());
  while (leaf->getNumChildren() > 0) 
  {
    leaf = leaf->getChild(0);
  }
  Genome* parent = leaf->getParent();  

  assert(firstTop <= lastTop);
  assert(lastTop < leaf->getNumTopSegments());
  assert(lastTop < parent->getNumBottomSegments());

  hal_index_t span = 1 + lastTop - firstTop;
  for (hal_index_t i = 0; i < span; ++i)
  {
    TopSegmentIteratorPtr top = leaf->getTopSegmentIterator(firstTop + i);
    BottomSegmentIteratorPtr bot =
       parent->getBottomSegmentIterator(lastTop - i);
    top->setParentIndex(NULL_INDEX);
    bot->setChildIndex(0, NULL_INDEX);
  }
}

void writeParseInfo(Genome* genome)
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
