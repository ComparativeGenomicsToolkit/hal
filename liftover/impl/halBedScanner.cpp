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

BedScanner::BedScanner() : _bedStream(NULL)
{

}

BedScanner::~BedScanner()
{

}

void BedScanner::scan(const string& bedPath, int bedVersion)
{
  assert(_bedStream == NULL);
  _bedStream = new ifstream(bedPath.c_str());
  
  try {
    scan(_bedStream, bedVersion);
  }
  catch(hal_exception e)
  {
    delete _bedStream;
    _bedStream = NULL;
    stringstream ss;
    ss << e.what() << " in file " << bedPath;
    throw hal_exception(ss.str());
  }  

  delete _bedStream;
  _bedStream = NULL;
}

void BedScanner::scan(istream* is, int bedVersion)
{
  _bedStream = is;
  if (bedVersion == -1)
  {
    bedVersion = getBedVersion(is);
  }

  if (_bedStream->bad())
  {
    throw hal_exception("Error reading bed input stream");
  }
  string lineBuffer;
  while (_bedStream->good())
  {
    _bedLine.read(*_bedStream, bedVersion, lineBuffer);
    visitLine();
  }
  visitEOF();
  _bedStream = NULL;
}

int BedScanner::getBedVersion(istream* bedStream)
{
  if (bedStream->bad())
  {
    throw hal_exception("Error reading bed input stream");
  }
  string lineBuffer;
  BedLine bedLine;
  int version = 12;
  streampos pos = bedStream->tellg();
  for (; version > 3; --version)
  {
    try
    {
      bedStream->seekg(pos);
      bedLine.read(*bedStream, version, lineBuffer);
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
  bedStream->seekg(pos);
  return version;
}

void BedScanner::visitLine()
{
}

void BedScanner::visitEOF()
{
}
