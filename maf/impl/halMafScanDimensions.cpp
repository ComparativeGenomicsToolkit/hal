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

  updateArrayIndices();
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
    pair<hal_size_t, hal_size_t> startIndex(0, (hal_size_t)-1);
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
      startIndex.first = 0;
      rec->_startMap.insert(startIndex);
    }
    rec->_length = row._srcLength;
    
    if (row._length > 0)
    {
      // add the begnning of the line as a segment start position
      // also add the last + 1 segments as a start position if in range
      hal_size_t start = row._startPosition;
      if (row._strand == '-')
      {
        start = row._srcLength - 1 - row._startPosition;
      }
      startIndex.first = start;
      rec->_startMap.insert(startIndex);
      hal_size_t end = start + row._length;
      if (end < row._srcLength)
      {
        startIndex.first = end;
        rec->_startMap.insert(startIndex);
      }

      size_t numGaps = 0;

      for (size_t j = 0; j < length; ++j)
      {
        if (row._line[j] == '-')
        {
          ++numGaps;
        }
        // valid segmentation between j-1 and j:
        // we add the start coordinate of the segment beginning at 
        // j in forward segment coordinates.  
        else if (_mask[j] == true && (row._line[j] != '-'))
        {
          startIndex.first = start + j - numGaps;
          rec->_startMap.insert(startIndex);
        }
      }
    }
  }
}

void MafScanDimensions::updateArrayIndices()
{
  for (DimMap::iterator i = _dimMap.begin(); i != _dimMap.end(); ++i)
  {
    assert(i->second != NULL);
    // compute the array index for each start position in the map
    hal_size_t arrayIndex = 0;
    StartMap& startMap = i->second->_startMap;
    for (StartMap::iterator j = startMap.begin(); j != startMap.end(); ++j)
    {
      j->second = arrayIndex++;
      assert(j->first >= j->second);
    }
  }
}
