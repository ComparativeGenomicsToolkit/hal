/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <deque>
#include "halRandomData.h"


using namespace std;
using namespace hal;

static inline hal_size_t uniformInt(mt19937& rng, hal_size_t min, hal_size_t max)
{
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

static inline double uniformDbl(mt19937& rng, double min, double max)
{
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

static inline bool exponEvent(mt19937& rng, double mu)
{
    std::uniform_real_distribution<double> dist;
    return dist(rng) <= (1.0 - exp(-mu));
}

static inline char randDNA(mt19937& rng)
{
    hal_size_t i = uniformInt(rng, 0, 3);
  switch (i)
  {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  default: return 'T';
  }
}

static inline void mutateString(mt19937& rng, string& buffer, double branchLength)
{
  for (size_t i = 0; i < buffer.length(); ++i)
  {
    if (exponEvent(rng, branchLength))
    {
      buffer[i] = randDNA(rng);
    }
  }
}

void hal::createRandomAlignment(mt19937& rng,
                                Alignment* emptyAlignment,
                                double meanDegree,
                                double maxBranchLength,
                                hal_size_t maxGenomes,
                                hal_size_t minSegmentLength,
                                hal_size_t maxSegmentLength,
                                hal_size_t minSegments,
                                hal_size_t maxSegments)
{
  createRandomTree(rng, emptyAlignment,
                   meanDegree,
                   maxBranchLength,
                   maxGenomes);
  
  createRandomDimensions(rng, emptyAlignment,
                         minSegmentLength,
                         maxSegmentLength,
                         minSegments,
                         maxSegments);

  deque<string> bfQueue;
  bfQueue.push_front(emptyAlignment->getRootName());

  while (bfQueue.empty() == false)
  {
    Genome* genome = emptyAlignment->openGenome(bfQueue.back());
    bfQueue.pop_back();

    createRandomGenome(rng, emptyAlignment, genome);
  
    vector<string> childNames = emptyAlignment->getChildNames(genome->getName());
    for (size_t i = 0; i < childNames.size(); ++i)
    {
      bfQueue.push_front(childNames[i]);
    }

    Genome* parent = genome->getParent();
    if (parent != NULL)
    {
      emptyAlignment->closeGenome(parent);
    }
    emptyAlignment->closeGenome(genome);
  }
}
                           

void hal::createRandomTree(mt19937& rng,
                           Alignment* emptyAlignment,
                           double meanDegree,
                           double maxBranchLength,
                           hal_size_t maxGenomes)
{
  assert(emptyAlignment->getNumGenomes() == 0);
  
  emptyAlignment->addRootGenome("Genome_0");
  
  deque<string> bfQueue;
  bfQueue.push_front(emptyAlignment->getRootName());
  size_t genomeCount = 1;
  
  while (bfQueue.empty() == false)
  {
    Genome* genome = emptyAlignment->openGenome(bfQueue.back());
    bfQueue.pop_back();
    hal_size_t numChildren = (hal_size_t)(uniformDbl(rng, 0.0, 2.0 * meanDegree) + 0.5);
    if (genomeCount + numChildren >= maxGenomes)
    {
      numChildren = maxGenomes - genomeCount;
    }

    for (hal_size_t i = 0; i < numChildren; ++i)
    {
        string childName = "Genome_" + std::to_string(genomeCount++);
      emptyAlignment->addLeafGenome(childName,
                                    genome->getName(),
                                    uniformDbl(rng, 1e-5, maxBranchLength));
      bfQueue.push_front(childName);
    }
  }
}

void hal::createRandomDimensions(mt19937& rng,
                                 Alignment* alignment,
                                 hal_size_t minSegmentLength,
                                 hal_size_t maxSegmentLength,
                                 hal_size_t minSegments,
                                 hal_size_t maxSegments)
{
  deque<string> bfQueue;
  bfQueue.push_front(alignment->getRootName());

  while (bfQueue.empty() == false)
  {
    Genome* genome = alignment->openGenome(bfQueue.back());
    assert(genome != NULL);
    bfQueue.pop_back();
    
    Genome* parent = genome->getParent();
    hal_size_t botSegSize = uniformInt(rng, minSegmentLength, maxSegmentLength);
    hal_size_t numBottomSegments = uniformInt(rng, minSegments, maxSegments);
    hal_size_t length = numBottomSegments * botSegSize;
    hal_size_t topSegSize = 0;
    hal_size_t numTopSegments = 0;
    if (parent != NULL)
    {
      BottomSegmentIteratorPtr it = parent->getBottomSegmentIterator();
      const BottomSegment* bseg = it->getBottomSegment();
      topSegSize = bseg->getLength();
      numTopSegments = length / topSegSize;
      if (length % topSegSize != 0)
      {
        ++numTopSegments;
      }
    }
    vector<string> childNames = alignment->getChildNames(genome->getName());
    if (childNames.empty())
    {
      numBottomSegments = 0;
    }
    if (numBottomSegments == 0 && numTopSegments == 0)
    {
      length = 0;
    }

    vector<Sequence::Info> dimensions;
    dimensions.push_back(Sequence::Info(genome->getName() + "_seq", 
                                        length, numTopSegments,
                                        numBottomSegments));
    
    genome->setDimensions(dimensions);

    hal_size_t numChildren = genome->getNumChildren();
    BottomSegmentIteratorPtr botIt = genome->getBottomSegmentIterator();
    for (size_t i = 0; i < genome->getNumBottomSegments(); ++i)
    {
      BottomSegment* bottomSegment = botIt->getBottomSegment();
      bottomSegment->setCoordinates(i * botSegSize, botSegSize);
      for (size_t j = 0; j < numChildren; ++j)
      {
        bottomSegment->setChildIndex(j, NULL_INDEX);
      }
      if (numTopSegments > 0)
      {
        bottomSegment->setTopParseIndex((i * botSegSize) / topSegSize);
      }
      else
      {
        bottomSegment->setTopParseIndex(NULL_INDEX);
      }
      botIt->toRight();
    }

    TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
    for (size_t i = 0; i < genome->getNumTopSegments(); ++i)
    {
      TopSegment* topSegment = topIt->getTopSegment();
      topSegment->setParentIndex(NULL_INDEX);
      hal_size_t tsLength = 0;
      if (i < genome->getNumTopSegments() - 1 || length % topSegSize == 0)
      {
        tsLength = topSegSize;
      }
      else
      {
        tsLength = length % topSegSize;
      }     
      topSegment->setCoordinates(i * topSegSize, tsLength);
      if (numBottomSegments > 0)
      {
        topSegment->setBottomParseIndex((i * topSegSize) / botSegSize);
      }
      else
      {
        topSegment->setBottomParseIndex(NULL_INDEX);
      }
      topIt->toRight();
    }

    for (size_t i = 0; i < childNames.size(); ++i)
    {
      bfQueue.push_front(childNames[i]);
    }
  }
}

void hal::createRandomGenome(mt19937& rng, Alignment* alignment, Genome* genome)
{
  Genome* parent = genome->getParent();
  set<pair<hal_index_t, hal_index_t> > edgeSet;
  if (parent == NULL)
  {
    DNAIteratorPtr dnaIt = genome->getDNAIterator();
    hal_size_t length = genome->getSequenceLength();
    for (hal_size_t i = 0; i < length; ++i)
    {
      dnaIt->setChar(randDNA(rng));
      dnaIt->toRight();
    }
  }
  else
  {
    vector<string> parentChildNames = 
       alignment->getChildNames(parent->getName());
    hal_size_t indexInParent = parentChildNames.size() ;
    for (hal_size_t i = 0; i < parentChildNames.size(); ++i)
    {
      if (parentChildNames[i] == genome->getName())
      {
        indexInParent = i;
      }
    }
    assert(indexInParent < parentChildNames.size());
    double branchLength = alignment->getBranchLength(parent->getName(),
                                                     genome->getName());
    
    TopSegmentIteratorPtr topIter = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIter = parent->getBottomSegmentIterator();
    hal_size_t numTopSegs = genome->getNumTopSegments();
    for (hal_size_t i = 0; i < numTopSegs; ++i)
    {
        createRandomSegment(rng, genome, indexInParent, 
                          edgeSet, topIter, botIter, branchLength);
      topIter->toRight();
    }
  }
}

void hal::createRandomSegment(mt19937& rng,
                              Genome* genome, 
                              hal_size_t indexInParent,
                              set<pair<hal_index_t, hal_index_t> >& edgeSet, 
                              TopSegmentIteratorPtr topIter, 
                              BottomSegmentIteratorPtr botIter,
                              double branchLength)
{
  Genome* parent = genome->getParent();
  hal_size_t numTopSegs = genome->getNumTopSegments();
  hal_size_t numBotSegs = parent ? parent->getNumBottomSegments() : 0;
  string buffer; 
  TopSegment* topSegment = topIter->getTopSegment();
  
  // case 1: parent index same as child index
  hal_index_t parentIdx = topSegment->getArrayIndex();

  // case 2: random parent index (trasposition/duplication)
  if (parentIdx >= (hal_index_t)numBotSegs || exponEvent(rng, branchLength) == true)
  {
      parentIdx = uniformInt(rng, 0, numBotSegs - 1);
  }

  // case 3: null parent index (insertion)
  else if (exponEvent(rng, branchLength) == true &&
           exponEvent(rng, branchLength) == true)
  {
    parentIdx = NULL_INDEX;
  }

  // case 4: don't know the size of the last segment so don't bother
  if (parentIdx == (hal_index_t)numBotSegs - 1 || 
      topSegment->getArrayIndex() == (hal_index_t)numTopSegs - 1)
  {
    parentIdx = NULL_INDEX;
  }

  //
  if (genome->getParent() == NULL)
  {
    parentIdx = NULL_INDEX;
  }

  topSegment->setParentIndex(parentIdx);
  topSegment->setParentReversed(false);
  topSegment->setNextParalogyIndex(NULL_INDEX);

  if (parentIdx == NULL_INDEX)
  {
    buffer.resize(topSegment->getLength());
    for (size_t j = 0; j < buffer.length(); ++j)
    {
      buffer[j] = randDNA(rng);
    }
  }
  
  else
  {    
      bool reversed = exponEvent(rng, branchLength);
    topSegment->setParentReversed(reversed);
    botIter->toParent(topIter);
    botIter->getString(buffer);

    mutateString(rng, buffer, branchLength);
    BottomSegment* botSegment = botIter->getBottomSegment();
    assert(botSegment->getArrayIndex() == parentIdx);
    assert(botSegment->getLength() == topSegment->getLength());
    botSegment->setChildIndex(indexInParent, topSegment->getArrayIndex());
    botSegment->setChildReversed(indexInParent, reversed);

    set<pair<hal_index_t, hal_index_t> >::iterator setIt;
    pair<hal_index_t, hal_index_t> key(botSegment->getArrayIndex(), 0);
    setIt = edgeSet.lower_bound(key);
    set<pair<hal_index_t, hal_index_t> >::iterator firstIt = setIt;
    set<pair<hal_index_t, hal_index_t> >::iterator prevIt = setIt;
    while (setIt != edgeSet.end() && setIt->first == botSegment->getArrayIndex())
    {
      prevIt = setIt++;
    }

    if (prevIt != edgeSet.end() && prevIt->first == botSegment->getArrayIndex())
    {
      assert (prevIt->second < topSegment->getArrayIndex());
      assert (prevIt->second != NULL_INDEX);
      TopSegmentIteratorPtr paralogousIt = 
         genome->getTopSegmentIterator(prevIt->second);
      TopSegment* paralogousSegment = paralogousIt->getTopSegment();
      paralogousSegment->setNextParalogyIndex(topSegment->getArrayIndex());
      topSegment->setNextParalogyIndex(firstIt->second);
    }

    edgeSet.insert(pair<hal_index_t, hal_index_t>(
                     botSegment->getArrayIndex(),
                     topSegment->getArrayIndex()));
  }

  genome->setSubString(buffer, topSegment->getStartPosition(),
                       topSegment->getLength());                           
}                           

