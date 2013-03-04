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
  assert(sizeof(ArrayInfo) == sizeof(hal_size_t));
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
  size_t dotPos = row._sequenceName.find('.');
  if (dotPos == string::npos || dotPos == 0 || 
      dotPos == row._sequenceName.length() - 1)
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
    pair<hal_size_t, ArrayInfo> startIndex;
    startIndex.first = 0;
    startIndex.second._index = 0; 
    startIndex.second._count = 1;
    startIndex.second._written = 0;
    startIndex.second._empty = 0;
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
      startIndex.second._empty = 1;
      rec->_startMap.insert(startIndex);
      startIndex.second._empty = 0;
      rec->_numSegments = 0;
    }
    rec->_length = row._srcLength;
    
    if (row._length > 0)
    {
      // add the begnning of the line as a segment start position
      // also add the last + 1 segments as a start position if in range
      hal_size_t start = row._startPosition;
      hal_size_t end = start + row._length;
      if (row._strand == '-')
      {
        start = row._srcLength - 1 - (row._startPosition + row._length - 1);
        end = row._srcLength - row._startPosition;
      }
      startIndex.first = start;
      pair<StartMap::iterator, bool> smResult = 
         rec->_startMap.insert(startIndex);
      StartMap::iterator smIt = smResult.first;
      bool bad = false;

      // check for duplication / inconsistency:
      // 1) new interval lands on start position of existing interval
      // existing is unchanged but we don't do anything else. 
      if (smResult.second == false && smIt->second._empty == 0)
      {
        bad = true;
      }

      // 2) new interval overlaps with existing interval
      // set count to 0 if new, ignore otherwise
      StartMap::iterator next = smIt;
      ++next;
      while (next != rec->_startMap.end() && !bad)
      {
        if (next->second._count > 0)
        {
          if (smIt->first + row._length > next->first)
          {
            bad = true;
          }
          else
          {
            break;
          }
        }
        ++next;
      }

      // 3) new interval overlaps a previous interval that is not 
      // empty
      if (!bad && smResult.second == true && smIt != rec->_startMap.begin())
      {
        StartMap::iterator prev = smIt;
        --prev;
        while (!bad)
        {
          if (prev->second._count > 0)
          {
            if (prev->second._empty == 0)
            {
              bad = true;
            }
            else
            {
              break;
            }
          }
          if (prev == rec->_startMap.begin())
          {
            break;
          }
          --prev;
        }
      }

      if (bad == true)
      {
        if (smResult.second == true)
        {
          rec->_startMap.erase(smIt);
        }
        rec->_badPosSet.insert(FilePosition(_mafFile.tellg(), i));
      }
      else
      {
        smIt->second._empty = 0;
        assert(smIt->second._count == 1);
        if (end < row._srcLength)
        {
          startIndex.first = end;
          startIndex.second._empty = 1;
          startIndex.second._count = 1;
          pair<StartMap::iterator, bool> endResult = 
             rec->_startMap.insert(startIndex);              
        }

        size_t numGaps = 0;
        hal_size_t interiorStart;

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
            interiorStart = start + j - numGaps;
            if (interiorStart > start)
            {
              ++smIt->second._count;
            }
          }
        }
      }
    }
  }
}

void MafScanDimensions::updateArrayIndices()
{
  string curName, prevName;
  hal_size_t arrayIndex = 0;
  hal_size_t prevIndex = 0;
  for (DimMap::iterator i = _dimMap.begin(); i != _dimMap.end(); ++i)
  {
    assert(i->second != NULL);
    curName = genomeName(i->first);
    assert(!curName.empty());
    if (curName != prevName)
    {
      arrayIndex = 0;
      prevIndex = 0;
    }
    // compute the array index for each start position in the map
    StartMap& startMap = i->second->_startMap;
    for (StartMap::iterator j = startMap.begin(); j != startMap.end(); ++j)
    {
      assert(j->second._count > 0);
      j->second._index = arrayIndex;
      arrayIndex += j->second._count;
    }    
    i->second->_numSegments = arrayIndex - prevIndex;
    prevName = curName;
    prevIndex = arrayIndex;
  }
}
