/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hal.h"
#include "defaultMappedSegment.h"

using namespace std;
using namespace hal;

DefaultMappedSegment::DefaultMappedSegment(
  DefaultSegmentIteratorConstPtr source,
  DefaultSegmentIteratorConstPtr target)
 :
  DefaultSegmentIterator(target->getStartOffset(), target->getEndOffset(), 
                         target->getReversed()),
  _segment(source->getSegment())
{
  assert(source->getLength() == target->getLength());
}

DefaultMappedSegment::~DefaultMappedSegment()
{

}
   
SegmentPtr DefaultMappedSegment::getSegment()
{
  return _segment;
}

SegmentConstPtr DefaultMappedSegment::getSegment() const
{
  return _segment;
}

hal_size_t DefaultMappedSegment::getNumSegmentsInGenome() const
{
  const Genome* genome = _segment->getGenome();
  return _segment->isTop() ? genome->getNumTopSegments() : 
     genome->getNumBottomSegments();
}

//////////////////////////////////////////////////////////////////////////////
// MAPPED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
SlicedSegmentConstPtr DefaultMappedSegment::getSource() const
{
  return _source;
}

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

hal_size_t DefaultMappedSegment::map(const DefaultSegmentIterator* source,
                                     vector<MappedSegmentConstPtr>& results,
                                     const Genome* tgtGenome,
                                     const set<const Genome*>* genomesOnPath,
                                     bool doDupes)
{
  assert(source != NULL);
  return 0;
}

bool DefaultMappedSegment::mapUp(vector<MappedSegmentConstPtr>& results)
{
  /* if (isTop() == true)
  {
    const Genome* parent = getGenome()->getParent();
    assert(parent != NULL);
    BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
    
    
    toParent();
  }
  else
  {
    hal_index_t rightCutoff = getEndPosition();
    toParseUp();
    
    
    
    }*/ return false;
}
