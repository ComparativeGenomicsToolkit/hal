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

#include "halMafScanner.h"

using namespace std;
using namespace hal;


MafScanner::MafScanner()
{

}

MafScanner::~MafScanner()
{

}

void MafScanner::scan(const string& mafFilePath, const set<string>& targets)
{
  _targets = targets;
  _mafFile.open(mafFilePath.c_str());
  _numBlocks = 0;

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
      if (_rows > 0)
      {
        updateMask();
        aLine();
        ++_numBlocks;
        nextLine();
      }
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
      if (_mafFile.bad() || _mafFile.fail())
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

      if (_targets.size() > 1 && // (will always include reference) 
          _targets.find(genomeName(row._sequenceName)) == _targets.end())
      {
        // genome not in targets, pretend like it never happened. 
        --_rows;
      }
      else
      {
        sLine();
      }
    }
    else
    {
      nextLine();
    }
  }
  if (_rows > 0)
  {
    updateMask();
    ++_numBlocks;
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
    _mask.resize(length);
    fill(_mask.begin(), _mask.end(), false);

    // scan left to right
    for (size_t i = 1; i < length; ++i)
    {
      // scan up to down
      for (size_t j = 0; j < _rows && _mask[i] == false; ++j)
      {
        // beginning of gap run. add position of first gap to mask
        if (_block[j]._line[i] == '-' && _block[j]._line[i-1] != '-')
        {
          _mask[i] = true;
        }
        // end of gap run. add position of first non gap to mask
        else if (_block[j]._line[i] != '-' && _block[j]._line[i-1] == '-')
        {
          _mask[i] = true;
        }
      }
    }
  }
}

string MafScanner::genomeName(const string& fullName)
{
  size_t dotPos = fullName.find('.');
  assert(dotPos != string::npos && dotPos > 0 && 
         dotPos < fullName.length() - 1);
  return fullName.substr(0, dotPos);
}

string MafScanner::sequenceName(const string& fullName)
{
  size_t dotPos = fullName.find('.');
  assert(dotPos != string::npos && dotPos > 0 && 
         dotPos < fullName.length() - 1);
  return fullName.substr(dotPos + 1);
}

