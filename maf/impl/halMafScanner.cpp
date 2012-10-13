/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "halMafScanner.h"

using namespace std;
using namespace hal;


MafScanner::MafScanner()
{

}

MafScanner::~MafScanner()
{

}

void MafScanner::scan(const std::string& mafFilePath)
{
  _mafFile.open(mafFilePath.c_str());

  if (!_mafFile)
  {
    throw runtime_error("error opening path: " + mafFilePath);
  }
  
  _rows = 0;
  _block.clear();
  string buffer;
  while (!_mafFile.eof() && _mafFile.good())
  {
    buffer.clear();
    _mafFile >> buffer;
    if (buffer == "a")
    {
      updateMask();
      aLine();
      nextLine();
      _rows = 0;
    }
    else if (buffer == "s")
    {
      ++_rows;
      if (_rows > _block.size())
      {
        _block.resize(_rows);
      }
      Row& row = _block[_rows - 1];
      _mafFile >> row._sequenceName >> row._startPosition >> row._length 
               >> row._strand >> row._srcLength >> row._line;
      if (!_mafFile.good())
      {
        throw hal_exception("error parsing sequence " + row._sequenceName);
      }
      if (_rows > 1 && row._line.length() != _block[_rows - 2]._line.length())
      {
        stringstream ss;
        ss << "two lines in same block have different lengths: " 
           << row._sequenceName << " " << row._startPosition << " and "
           << _block[_rows - 2]._sequenceName << " " 
           << _block[_rows - 2]._startPosition;
        throw hal_exception(ss.str());
      }
      sLine();
    }
    else
    {
      nextLine();
    }
  }
  end();  
  _mafFile.close();
}

void MafScanner::nextLine()
{
  while (!_mafFile.eof() && !_mafFile.bad() && _mafFile.peek() != '\n')
  {
    _mafFile.get();
  }
}

// the mask stores a bit for every column where a gap begins in any row
// a mask at position i implies segmentation [0-i-1][i-n].  ie the cut is 
// on the left. 
void MafScanner::updateMask()
{
  if (_rows > 0)
  {
    size_t length = _block[0]._line.length();
    _mask.resize(length, false);
    
    for (size_t i = 1; i < length; ++i)
    {
      for (size_t j = 0; j < _rows && _mask[j] == true; ++j)
      {
        // beginning of gap run. add position of first gap to mask
        if (_block[j]._line[i] == '-' && _block[j]._line[i-1] != '-')
        {
          _mask[i] = true;
        }
        // end of gap run. add position of first non gap to mask
        else if (_block[j]._line[j] != '-' && _block[j]._line[i-1] == '-')
        {
          _mask[i] = true;
        }
      }
    }
  }
}

