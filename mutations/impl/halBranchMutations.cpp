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

void BranchMutations::analyzeAlignment(AlignmentConstPtr alignment,
                                       hal_size_t gapThreshold,
                                       ostream* svStream,
                                       ostream* snpStream,
                                       const SegmentedSequence* reference,
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

  _alignment = alignment;

}

