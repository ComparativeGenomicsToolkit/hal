/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halDetailedMutations.h"

using namespace std;
using namespace hal;

DetailedMutations::DetailedMutations()
{

}

DetailedMutations::~DetailedMutations()
{

}

void DetailedMutations::analyzeAlignment(AlignmentConstPtr alignment,
                                         hal_size_t gapThreshold,
                                         ostream* snpStream,
                                         ostream* svStream,
                                         const SegmentedSequence* reference,
                                         hal_index_t startPosition,
                                         hal_size_t length,
                                         const set<const Genome*>* targets)
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
  hal_index_t lastPosition = startPosition + (hal_index_t)(length - 1);

  _alignment = alignment;

  ColumnIteratorConstPtr colIt = reference->getColumnIterator(targets,
                                                              gapThreshold, 
                                                              startPosition,
                                                              lastPosition);
  
  while (colIt->lastColumn() == false)
  {
    colIt->toRight();
    
    // erase empty entries from the column.  helps when there are 
    // millions of sequences (ie from fastas with lots of scaffolds)
    // if (numBlocks++ % 1000 == 0)
    // {
    //  colIt->defragment();
    // }
  }

}

