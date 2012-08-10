/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <deque>
#include "defaultRearrangement.h"
#include "hal.h"

using namespace std;
using namespace hal;

// maximum size a simple indel can be to be considered a gap (and not
// a rearrangement)
const hal_size_t DefaultRearrangement::DefaultGapThreshold = 10;

DefaultRearrangement::DefaultRearrangement() :
  _gapThreshold(DefaultGapThreshold),
  _id(Invalid),
  _length(0),
  _numGaps(0),
  _numGapBases(0),
  _dupDegree(0)
{

}

DefaultRearrangement::~DefaultRearrangement()
{

}
   
// Rearrangement Interface Methods
DefaultRearrangement::
ID DefaultRearrangement::getID() const
{
  return _id;
}

hal_size_t DefaultRearrangement::getLength() const
{
  return _length;
}

hal_size_t DefaultRearrangement::getNumContainedGaps() const
{
  return _numGaps;
}

hal_size_t DefaultRearrangement::getNumContainedGapBases() const
{
  return _numGapBases;
}

hal_size_t DefaultRearrangement::getDuplicationDegree() const
{
  return _dupDegree;
}

TopSegmentIteratorConstPtr DefaultRearrangement::getLeftBreakpoint() const
{
  return _left;
}

TopSegmentIteratorConstPtr DefaultRearrangement::getRightBreakpoint() const
{
  return _right;
}

bool DefaultRearrangement::identifyFromLeftBreakpoint(
  TopSegmentIteratorConstPtr topSegment)
{
  return false;
}

bool DefaultRearrangement::identifyNext()
{
  return false;
}

hal_size_t DefaultRearrangement::getGapLengthThreshold() const
{
  return _gapThreshold;
}

void DefaultRearrangement::setGapLengthThreshold(hal_size_t threshold)
{
  _gapThreshold = threshold;
}
