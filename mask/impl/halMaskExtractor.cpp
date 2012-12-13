/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include "halMaskExtractor.h"

using namespace std;
using namespace hal;


MaskExtractor::MaskExtractor()
{

}

MaskExtractor::~MaskExtractor()
{

}

void MaskExtractor::extract(AlignmentConstPtr alignment, ostream* bedStream, 
                            hal_size_t extend, double extendPct)
{
  _alignment = alignment;
  _bedStream = bedStream;
  _extend = extend;
  _extendPct = extendPct;

}
