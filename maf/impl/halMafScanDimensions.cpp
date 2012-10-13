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
#include "halMafScanDimensions.h"

using namespace std;
using namespace hal;


MafScanDimensions::MafScanDimensions() : MafScanner()
{

}

MafScanDimensions::~MafScanDimensions()
{
  for (DimMap::iterator i = _dimMap.begin(); i != _dimMap.end(); ++i)
  {
    delete i->second;
  }
}

void MafScanDimensions::scan(const string& mafPath, const set<string>& targets)
{
  for (DimMap::iterator i = _dimMap.begin(); i != _dimMap.end(); ++i)
  {
    delete i->second;
  }
  _dimMap.clear();

  MafScanner::scan(mafPath, targets);

  updateDimensionsGlobal();
}

const MafScanDimensions::DimMap& MafScanDimensions::getDimensions() const
{
  return _dimMap;
}

void MafScanDimensions::aLine()
{
  assert(_rows <= _block.size());
  if (_rows > 0)
  {
    updateDimensionsFromBlock();
  }
}

void MafScanDimensions::sLine()
{
  Row& row = _block[_rows - 1];

  // this is the first pass.  so we do a quick sanity check
  if (row._sequenceName.find('.') == string::npos)
  {
    stringstream ss;
    ss << "illegal sequence name found: " << row._sequenceName << ".  Sequence "
       "names must be in genomeName.sequenceName format.";
    throw hal_exception(ss.str());
  }
  
  size_t numGaps = 0;
  for (size_t i = 0; i < row._line.length(); ++i)
  {
    if (row._line[i] == '-')
    {
      ++numGaps;
    }
    else
    {
      if (!isNucleotide(row._line[i]))
      {
        stringstream ss;
        ss << "problem reading line for sequence " << row._sequenceName << ": "
           " non-gap non-nucleotide character " << row._line[i] << " found at "
           << i << "th position";
        throw hal_exception(ss.str());
      }
    }
  }

  if (row._line.length() - numGaps != row._length)
  {
    stringstream ss;
    ss << "problem reading line for sequence " << row._sequenceName << ": "
           " length field was " << row._length << " but line contains "
           << row._line.length() - numGaps << " bases";
    throw hal_exception(ss.str());
  }

  if (row._startPosition + row._length > row._srcLength)
  {
    stringstream ss;
    ss << "problem reading line for sequence " << row._sequenceName << ": "
           " sequence length is" << row._srcLength << " but line starts at "
       << row._startPosition << " and contains " << row._length << " bases";
    throw hal_exception(ss.str());
  }
}

void MafScanDimensions::end()
{
  assert(_rows <= _block.size());
  if (_rows > 0)
  {
    updateDimensionsFromBlock();
  }
}

void MafScanDimensions::updateDimensionsFromBlock()
{
  assert(_rows > 0 && !_block.empty());
  size_t length = _block[0]._line.length();
  for (size_t i = 0; i < _rows; ++i)
  {
    Row& row = _block[i];
    pair<string, Record*> newRec(row._sequenceName, NULL);
    pair<DimMap::iterator, bool> result = _dimMap.insert(newRec);
    Record*& rec = result.first->second;
    if (result.second == false && row._srcLength != rec->_length)
    {
      assert(rec != NULL);
      stringstream ss;
      ss << "conflicting length for sequence " << row._sequenceName << ": "
         << "was scanned once as " << row._srcLength << " then again as "
         << rec->_length;
      throw hal_exception(ss.str());
    }
    else if (result.second == true)
    {
      rec = new Record(); 
      rec->_segments = 0;
    }
    rec->_length = row._srcLength;
    
    if (row._length > 0)
    {
      // any non-empty block increases are segment count by at least one.
      ++rec->_segments;
      
      // add the interval to the list.  we need to store this to infer
      // segments for intervals that *are not* covered by blocks. 
      rec->_intervals.push_back(Interval());
      Interval& ival = rec->_intervals.back();
      if (row._strand == '+')
      {
        ival.first = row._startPosition;
        ival.second = row._startPosition + row._length - 1;
      }
      else
      {
        ival.second = row._srcLength - 1 - row._startPosition;
        ival.first = ival.second - row._length + 1; 
      }
      assert(ival.first <= ival.second);
    }
    for (size_t j = 1; j < length; ++j)
    {
      // valid segmentation between j-1 and j
      if (_mask[j] == true && (row._line[j-1] != '-' || row._line[j] != '-'))
      {
        ++rec->_segments;
      }
    }
  }
}

void MafScanDimensions::updateDimensionsGlobal()
{
  for (DimMap::iterator i = _dimMap.begin(); i != _dimMap.end(); ++i)
  {
    Record*& rec = i->second;
    sort(rec->_intervals.begin(), rec->_intervals.end());

    if (!rec->_intervals.empty())
    {
      Interval& first = rec->_intervals.at(0);
      if (first.first > 0)
      {
        // empty segment before first interval
        ++rec->_segments;
      }
      Interval& last = rec->_intervals.back();
      if (last.second < rec->_length - 1)
      {
        // empty segment after last interval
        ++rec->_segments;
      }
    }

    for (size_t j = 1; j < rec->_intervals.size(); ++j)
    {
      Interval& cur = rec->_intervals.at(j);
      Interval& prev = rec->_intervals.at(j-1);
      if (prev.second >= cur.first)
      {
        stringstream ss;
        ss << "overlapping blocks detected for " << i->first 
           << "(" << cur.first << "," << cur.second << ") and ("
           << prev.first << "," << prev.second << ")";
        throw hal_exception(ss.str());
      }
      if (cur.first - prev.second > 1)
      {
        // empty segment between prev and cur
        ++rec->_segments;
      }
    }
  }
}
