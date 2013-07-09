/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halWiggleLiftover.h"

using namespace std;
using namespace hal;

WiggleLiftover::WiggleLiftover()
{

}

WiggleLiftover::~WiggleLiftover()
{

}

void WiggleLiftover::convert(AlignmentConstPtr alignment,
                             const Genome* srcGenome,
                             istream* inputFile,
                             const Genome* tgtGenome,
                             ostream* outputFile,
                             bool traverseDupes,
                             bool unique)
{
  _alignment = alignment;
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _outStream = outputFile;
  _traverseDupes = traverseDupes;
  _unique = unique;

  scan(inputFile);
}

void WiggleLiftover::visitHeader()
{
  _srcSequence = _srcGenome->getSequence(_sequenceName);
}

void WiggleLiftover::visitLine()
{

}
                      
