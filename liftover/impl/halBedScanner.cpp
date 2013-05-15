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
  if (bedVersion == -1)
  {
    bedVersion = getBedVersion(bedPath);
  }

  _bedFile.open(bedPath.c_str());
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
  visitEOF();
  _bedFile.close();
}

int BedScanner::getBedVersion(const string& bedPath)
{
  ifstream bedFile(bedPath.c_str());
  if (bedFile.bad())
  {
    stringstream ss;
    ss << "Error opening bed file: " << bedPath;
    throw hal_exception(ss.str());
  }
  string lineBuffer;
  BedLine bedLine;
  int version = 12;
  for (; version > 3; --version)
  {
    try
    {
      bedLine.read(bedFile, version, lineBuffer);
      break;
    }
    catch(...)
    {
      if (version == 12)
      {
        version = 10;
      }
    }
  }
  
  bedFile.close();
  return version;
}

void BedScanner::visitLine()
{
}

void BedScanner::visitEOF()
{
}
