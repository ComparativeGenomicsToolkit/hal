/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halBranchMutations.h"

using namespace std;
using namespace hal;

BranchMutations::BranchMutations()
{

}

BranchMutations::~BranchMutations()
{

}

void BranchMutations::analyzeBranch(AlignmentConstPtr alignment,
                                    hal_size_t gapThreshold,
                                    ostream* refBedStream,
                                    ostream* delBedStream,
                                    ostream* snpBedStream,
                                    const Genome* reference,
                                    hal_index_t startPosition,
                                    hal_size_t length)
{
  assert(reference != NULL);
  if (startPosition >= (hal_index_t)reference->getSequenceLength() ||
      (hal_size_t)startPosition + length > reference->getSequenceLength())
  {
    throw hal_exception("Invalid range specified for convertGenome");
  }
  if (length == 0)
  {
    length = reference->getSequenceLength() - startPosition;
  }
  if (length == 0)
  {
    throw hal_exception("Cannot convert zero length sequence");
  }

  _reference = reference;
  _start = startPosition;
  _length = length;
  _alignment = alignment;
  _maxGap = gapThreshold;
  _refStream = refBedStream;
  _delStream = delBedStream;
  _snpStream = snpBedStream;

  hal_index_t end = startPosition + (hal_index_t)length - 1;

  TopSegmentIteratorConstPtr top  = reference->getTopSegmentIterator();
  top->toSite(startPosition);
  
  _rearrangement = reference->getRearrangement(top->getArrayIndex());
  _rearrangement->setGapLengthThreshold(_maxGap);
  
  do {
    switch (_rearrangement->getID())
    {
    case Rearrangement::Inversion:     
    case Rearrangement::Deletion:
    case Rearrangement::Transposition:
    case Rearrangement::Insertion:
    case Rearrangement::Duplication:
      writeRearrangement();
    default:
      break;
    }
    writeSubstitutions();
    writeGapInsertions();
  }
  while (_rearrangement->identifyNext() == true && 
         _rearrangement->getLeftBreakpoint()->getStartPosition() <= end);
  
  // kind of stupid, but we do a second pass to get the gapped deletions
  _rearrangement = reference->getRearrangement(top->getArrayIndex());   
  _rearrangement->setAtomic(true);
  do {
    switch (_rearrangement->getID())
    {
    case Rearrangement::Deletion:
      writeGapDeletion();
    default:
      break;
    }
  }
  while (_rearrangement->identifyNext() == true && 
         _rearrangement->getLeftBreakpoint()->getStartPosition() <= end);
}

void BranchMutations::writeRearrangement()
{

}

void BranchMutations::writeSubstitutions()
{

}

void BranchMutations::writeGapInsertions()
{

}

void BranchMutations::writeGapDeletion()
{

}
