/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "halMafScanReference.h"

using namespace std;
using namespace hal;


MafScanReference::MafScanReference() : MafScanner()
{
}

MafScanReference::~MafScanReference()
{
}

std::string MafScanReference::getRefName(const std::string& mafPath)
{
  MafScanner::scan(mafPath, set<string>());
  return _name;
}

void MafScanReference::aLine()
{
}

void MafScanReference::sLine()
{
  Row& row = _block[_rows - 1];
  // this is the first pass.  so we do a quick sanity check
  if (row._sequenceName.find('.') == string::npos || 
      row._sequenceName.find('.') == 0)
  {
    stringstream ss;
    ss << "illegal sequence name found: " << row._sequenceName << ".  Sequence "
       "names must be in genomeName.sequenceName format.";
    throw hal_exception(ss.str());
  }
  
  _name = genomeName(row._sequenceName);
  _mafFile.seekg(0, ios_base::end);
}

void MafScanReference::end()
{
}
