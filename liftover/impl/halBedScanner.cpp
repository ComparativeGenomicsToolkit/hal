/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>

#include "halBedScanner.h"

using namespace std;
using namespace hal;

BedScanner::BedScanner() 
{

}

BedScanner::~BedScanner()
{

}

void BedScanner::scan(const string& bedPath, int bedVersion)
{
  _bedFile.open(bedPath.c_str());
  if (bedVersion == -1)
  {
    bedVersion = getBedVersion();
  }
  if (_bedFile.bad())
  {
    stringstream ss;
    ss << "Error opening bed file: " << bedPath;
    throw hal_exception(ss.str());
  }
  string lineBuffer;
  while (_bedFile.good())
  {
    _bedLine.read(_bedFile, bedVersion, lineBuffer);
    visitLine();
  }
}

