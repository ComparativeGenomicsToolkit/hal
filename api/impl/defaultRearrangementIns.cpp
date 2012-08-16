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

bool 
DefaultRearrangement::scanInsertionCycle(TopSegmentIteratorConstPtr topSegment)
{
  return false;
}

bool DefaultRearrangement::middleInsertion()
{
  return false;
}

bool DefaultRearrangement::startInsertion()
{
  return false;
}

bool DefaultRearrangement::endInsertion()
{
  return false;
}
