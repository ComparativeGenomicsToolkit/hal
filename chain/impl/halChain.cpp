/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halChain.h"

using namespace std;
using namespace hal;

static void convertBlocks(TopSegmentIteratorConstPtr firstIt, 
                          Chain& outChain);

ostream& hal::operator<<(ostream& os, const ChainBlock& b)
{
  return os << b._size << ' ' << b._tGap << ' ' << b._qGap;
}

ostream& hal::operator<<(ostream& os, const Chain& c)
{
  // score hardcoded to 0 for now
  os << "chain 0 " 
     << c._tName << ' '
     << c._tSize << ' '
     << c._tStrand << ' '
     << c._tStart << ' '
     << c._tEnd << ' '
     << c._qName << ' '
     << c._qSize << ' '
     << c._qStrand << ' '
     << c._qStart << ' '
     << c._qEnd << ' '
     << c._id << '\n';
  
  for (size_t i = 0; i < c._blocks.size(); ++i)
  {
    os << c._blocks[i] << '\n';
  }
  
  return os;
}

void hal::gtIteratorToChain(GappedTopSegmentIteratorConstPtr top, 
                            Chain& outChain,
                            hal_offset_t startOffset, 
                            hal_offset_t endOffset)
{
  hal_size_t childIndex = top->getChildIndex();
  const Genome* qGenome = top->getGenome();
  const Genome* tGenome = qGenome->getParent();
  assert(tGenome != NULL);

  // should be an input parameter to avoid recreating each time.
  GappedBottomSegmentIteratorConstPtr bottom = 
     tGenome->getGappedBottomSegmentIterator(0, childIndex, 
                                             top->getGapThreshold(), 
                                             top->getAtomic());

  bottom->toParent(top);

  // set up global information
  const Sequence* qSequence = top->getSequence();
  const Sequence* tSequence = bottom->getSequence();
  
  outChain._tName = tSequence->getName();
  outChain._qName = qSequence->getName();
  outChain._tSize = tSequence->getSequenceLength();
  outChain._qSize = qSequence->getSequenceLength();

  // slice iterator coordinates according to given offsets
  outChain._qStart = top->getStartPosition() + startOffset;
  outChain._qEnd = top->getStartPosition() - endOffset;
  outChain._qStrand = '+';
  outChain._qStart -= qSequence->getStartPosition();
  outChain._qEnd -= qSequence->getStartPosition();
  assert(outChain._qEnd > outChain._qStart);

  outChain._tStart = bottom->getStartPosition() + startOffset;
  outChain._tEnd = bottom->getEndPosition() - endOffset;
  outChain._tStart -= tSequence->getStartPosition();
  outChain._tEnd -= tSequence->getStartPosition();

  if (top->getParentReversed() == true)
  {
    outChain._tStart = outChain._tSize - outChain._tStart;
    outChain._tEnd = outChain._tSize - outChain._tEnd;
    outChain._tStrand = '-';
  }
  else
  {
    outChain._tStrand = '+';
  }
  assert(outChain._tEnd > outChain._tStart);
  
  // convert blocks
  TopSegmentIteratorConstPtr firstIt = top->getLeft();
  convertBlocks(firstIt, outChain);
}

void convertBlocks(TopSegmentIteratorConstPtr firstIt, 
                   Chain& outChain)
{
  TopSegmentIteratorConstPtr topIt = firstIt->copy();
  const Sequence* sequence = topIt->getSequence();
  hal_index_t start = (hal_index_t)outChain._tStart -
     sequence->getStartPosition();
  hal_index_t end = (hal_index_t)outChain._tEnd - sequence->getStartPosition();
  BottomSegmentIteratorConstPtr botIt =
     topIt->getGenome()->getParent()->getBottomSegmentIterator();

  hal_index_t tPrev = NULL_INDEX;
  hal_index_t qPrev = NULL_INDEX;
  
  vector<ChainBlock>& blocks = outChain._blocks;
  blocks.resize(0);

  while (topIt->getEndPosition() >= start && 
         topIt->getStartPosition() < end)
  {
    if (topIt->getStartOffset() == 0 && topIt->getEndOffset() == 0)
    {
      hal_size_t startOffset = 0;
      hal_size_t endOffset = 0;

      if (start > topIt->getStartPosition() && start < topIt->getEndPosition())
      {
        startOffset = start - topIt->getStartPosition();
      }
      if (end-1 > topIt->getStartPosition() && end-1 < topIt->getEndPosition())
      {
        endOffset = topIt->getEndPosition() - end + 1;
      }
      topIt->slice(startOffset, endOffset);
    }
       
    if (topIt->hasParent() == true)
    {
      botIt->toParent(topIt);
      blocks.resize(blocks.size() + 1);
      blocks.back()._size = topIt->getLength();
      blocks.back()._tGap = 0;
      blocks.back()._qGap = 0;
 
      if (tPrev != NULL_INDEX)
      {
        blocks[blocks.size() - 2]._tGap = topIt->getStartPosition() - tPrev;
        blocks[blocks.size() - 2]._qGap = std::abs(botIt->getStartPosition() -
                                                   qPrev);
      }
      
      tPrev = topIt->getEndPosition();
      qPrev = botIt->getEndPosition();
    }
    topIt->toRight();
  }
}
  
