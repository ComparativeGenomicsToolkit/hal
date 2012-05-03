/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "defaultColumnIterator.h"
#include "hal.h"

using namespace std;
using namespace hal;

DefaultColumnIterator::DefaultColumnIterator(Genome* reference, 
                                             Genome* root,
                                             hal_index_t columnIndex)
{

}
   
DefaultColumnIterator::~DefaultColumnIterator()
{

}

void DefaultColumnIterator::toRight() const
{

}

bool DefaultColumnIterator::equals(ColumnIteratorConstPtr other) const
{
  return false;
}

const ColumnIterator::SegmentMap& DefaultColumnIterator::getSegmentMap() const
{
  return _segmentMap;
}

const Genome* ColumnIterator::getReferenceGenome() const 
{
  return NULL;
}

