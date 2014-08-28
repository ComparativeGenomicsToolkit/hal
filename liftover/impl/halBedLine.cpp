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

#include "halBedLine.h"

using namespace std;
using namespace hal;

BedLine::BedLine() : _start(NULL_INDEX), _end(NULL_INDEX), _strand('+'),
                     _version(NULL_INDEX), _srcStart(NULL_INDEX)
{

}

BedLine::~BedLine()
{

}

istream& BedLine::read(istream& is, int version, string& lineBuffer)
{
  _version = version;
  std::getline(is, lineBuffer);
  stringstream ss(lineBuffer);
  ss.imbue(is.getloc());
  ss >> _chrName;
  if (ss.bad() || ss.fail()) 
  {
    throw hal_exception("Error scanning BED chrom");
  }
  ss >> _start;
  if (ss.bad() || ss.fail()) 
  {
    throw hal_exception("Error scanning BED chromStart");
  }
  ss >> _end;
  if (ss.bad() || ss.fail())
  {
    throw hal_exception("Error scanning BED chromEnd");
  }
  if (_version > 3)
  {
    ss >> _name;
    if (ss.bad() || ss.fail())
    {
      throw hal_exception("Error scanning BED name");
    }
  }
  if (_version > 4)
  {
    ss >> _score;
    if (ss.bad() || ss.fail())
    {
      throw hal_exception("Error scanning BED score");
    }
  }
  if (_version > 5)
  {
    ss >> _strand;
    if (ss.bad() || ss.fail())
    {
      throw hal_exception("Error scanning BED strand");
    }
    if (_strand != '.' && _strand != '+' && _strand != '-')
    {
      throw hal_exception("Strand character must be + or - or .");
    }
  }
  if (_version > 6)
  {
    ss >> _thickStart;
    if (ss.bad() || ss.fail())
    {
      throw hal_exception("Error scanning BED thickStart");
    }
  }
  if (_version > 7)
  {
    ss >> _thickEnd;
    if (ss.bad() || ss.fail()) 
    {
      throw hal_exception("Error scanning BED thickEnd");
    }
  }
  if (_version > 8)
  {
    string rgb;
    ss >> rgb;    
    if (ss.bad() || ss.fail())
    {
      throw hal_exception("Error scanning BED itemRGB");
    }
    vector<string> rgbTokens = chopString(rgb, ",");
    if (rgbTokens.size() > 3 || rgbTokens.size() == 0)
    {
      throw hal_exception("Error parsing BED itemRGB");
    }
    stringstream rgbssr(rgbTokens[0]);
    rgbssr >> _itemR;
    if (rgbssr.bad())
    {
      throw hal_exception("Error parsing BED itemRGB");
    }
    _itemG = _itemR;
    _itemB = _itemR;
    if (rgbTokens.size() > 1)
    {
      stringstream rgbssg(rgbTokens[1]);
      rgbssg >> _itemG;
      if (rgbssg.bad())
      {
        throw hal_exception("Error parsing BED itemRGB");
      }
    }
    if (rgbTokens.size() == 3)
    {
      stringstream rgbssb(rgbTokens[2]);
      rgbssb >> _itemB;
      if (rgbssb.bad())
      {
        throw hal_exception("Error parsing BED itemRGB");
      }    
    }
  }
  if (_version > 9)
  {
    size_t numBlocks;
    ss >> numBlocks;
    if (ss.bad() || ss.fail()) 
    {
      throw hal_exception("Error scanning BED blockCount");
    }
    if (numBlocks > 0)
    {
      string blockSizes;
      ss >> blockSizes;
      if (ss.bad() || ss.fail())
      {
        throw hal_exception("Error scanning BED blockSizes");
      }
      string blockStarts;
      ss >> blockStarts;
      if (ss.bad() || ss.fail())
      {
        throw hal_exception("Error scanning BED blockStarts");
      }
      _blocks.resize(numBlocks);
      vector<string> sizeBuf = chopString(blockSizes, ",");
      if (sizeBuf.size() != numBlocks)
      {
        throw hal_exception("Error scanning BED blockSizes");
      }
      vector<string> startBuf = chopString(blockStarts, ",");
      if (startBuf.size() != numBlocks)
      {
        throw hal_exception("Error scanning BED blockStarts");
      }
      for (size_t i = 0; i < numBlocks; ++i)
      {
        BedBlock& block = _blocks[i];
        stringstream ss1(sizeBuf[i]);
        ss1 >> block._length;
        if (ss1.bad())
        {
          throw hal_exception("Error scanning BED blockSizes");
        }
        stringstream ss2(startBuf[i]);
        ss2 >> block._start;
        if (ss2.bad())
        {
          throw hal_exception("Error scanning BED blockStarts");
        }
        if (_start + block._start + block._length > _end)
        {
          throw hal_exception("Error BED block out of range");
        }
      }
    }
  }
  _extra.clear();
  while (ss.good())
  {
    string extraBuf;
    ss >> extraBuf;
    if (extraBuf.length() > 0)
    {
      _extra.push_back(extraBuf);
    }
  }
  return is;
}

ostream& BedLine::write(ostream& os, int version)
{
  if (version <= 0)
  {
    version = _version;
  }
  os << _chrName << '\t' << _start << '\t' << _end;
  
  if (version > 3)
  {
    os << '\t' << _name;
  }
  if (version > 4)
  {
    os << '\t' << _score;
  }
  if (version > 5)
  {
    os << '\t' << _strand;
  }
  if (version > 6)
  {
    os << '\t' << _thickStart;
  }
  if (version > 7)
  {
    os << '\t' << _thickEnd;
  }
  if (version > 8)
  {
    os << '\t' << _itemR << ',' << _itemG << ',' << _itemB;
  }
  if (version > 9)
  {
    os << '\t' << _blocks.size();
    for (size_t i = 0; i < _blocks.size(); ++i)
    {
      if (i == 0)
      {
        os << '\t';
      }
      else
      {
        os << ',';
      }
      os << _blocks[i]._length;
    }
    for (size_t i = 0; i < _blocks.size(); ++i)
    {
      if (i == 0)
      {
        os << '\t';
      }
      else
      {
        os << ',';
      }
      os << _blocks[i]._start;
    }
  }

  for (size_t i = 0; i < _extra.size(); ++i)
  {
    os << '\t' << _extra[i];
  }
  
  os << '\n';
  return os;
}

bool BedBlock::operator<(const BedBlock& other) const
{
  bool isLess = _start < other._start;
  return isLess;
}

bool BedLineLess::operator()(const BedLine& b1, const BedLine& b2) const 
{
  if (b1._chrName < b2._chrName)
  {
    return true;
  }
  else if (b1._chrName == b2._chrName)
  {
    if (b1._start < b2._start)
    {
      return true;
    }
    else if (b1._start == b2._start)
    {
      if (b1._end < b2._end)
      {
        return true;
      }
      else if (b1._end == b2._end)
      {
        return b1._strand == '+' && b2._strand != '+';
      }
    }
  }
  return false;
}

bool BedLineSrcLess::operator()(const BedLine& b1, const BedLine& b2) const 
{
  return b1._srcStart < b2._srcStart;
}


ostream& BedLine::writePSL(ostream& os, bool prefixWithName)
{
  assert(_psl.size() == 1);
  const PSLInfo& psl = _psl[0];
  assert(_blocks.size() == psl._qBlockStarts.size());
  assert(_blocks.size() > 0);
  assert(_srcStart >= (hal_index_t)psl._qChromOffset);
  if (validatePSL() == false)
  {
    throw hal_exception("Internal error: PSL does not validate");
  }

  if (prefixWithName == true)
  {
    os << _name << '\t';
  }
  os << psl._matches << '\t'
     << psl._misMatches << '\t'
     << psl._repMatches << '\t'
     << psl._nCount << '\t'
     << psl._qNumInsert << '\t'
     << psl._qBaseInsert << '\t'
     << psl._tNumInsert << '\t'
     << psl._tBaseInsert << '\t'
     << psl._qStrand << _strand << '\t'
     << psl._qSeqName << '\t'
     << psl._qSeqSize << '\t'
     << (_srcStart - psl._qChromOffset) << '\t'
     << (psl._qEnd - psl._qChromOffset) << '\t'
     << _chrName << '\t'
     << psl._tSeqSize << '\t'
     << _start << '\t'
     << _end << '\t'
     << _blocks.size() << '\t';

  for (size_t i = 0; i < _blocks.size(); ++i)
  {
    os << _blocks[i]._length << ',';
  }
  os << '\t';

  for (size_t i = 0; i < psl._qBlockStarts.size(); ++i)
  {
    assert(psl._qBlockStarts[i] >= (hal_index_t)psl._qChromOffset);
    hal_index_t start = psl._qBlockStarts[i] - psl._qChromOffset;
    if (psl._qStrand == '-')
    {
      start = psl._qSeqSize - start - _blocks[i]._length;
    }
    os << start << ',';
  }
  os << '\t';

  for (size_t i = 0; i < _blocks.size(); ++i)
  {
    hal_index_t start = _blocks[i]._start + _start;
    if (_strand == '-')
    {
      start = psl._tSeqSize - start - _blocks[i]._length;
    }
    os << start << ',';
  }

  os << '\n';
  return os;
}

bool BedLine::validatePSL() const
{
  if (_psl.size() != 1)
  {
    assert(false); return false;
  }
  const PSLInfo& psl = _psl[0];

  if (_blocks.size() < 1)
  {
    assert(false); return false;
  }
  if (_blocks.size() != psl._qBlockStarts.size())
  {
    assert(false); return false;
  }

  if (psl._qBlockStarts.size() != _blocks.size())
  {
    assert(false); return false;
  }

  hal_size_t totBlockLen = 0;
  for (size_t i = 0; i < _blocks.size(); ++i)
  {
    totBlockLen += _blocks[i]._length;
  }
  
  if (totBlockLen != psl._matches + psl._misMatches + psl._repMatches)
  {
    assert(false); return false;
  }

  if (totBlockLen + psl._qBaseInsert != psl._qEnd - _srcStart)
  {
    assert(false); return false;
  }
  
  if (totBlockLen + psl._tBaseInsert != (hal_size_t)_end - _start)
  {
    assert(false); return false;
  }

  if (_strand != '-')
  {
    if (_blocks[0]._start != 0)
    {
      assert(false); return false;
    }
    if (_blocks.back()._start + _blocks.back()._length + _start != _end)
    {
      assert(false); return false;
    }
  }
  else
  {
    if (_blocks.back()._start != 0)
    {
      assert(false); return false;
    }
    if (_blocks[0]._start + _blocks[0]._length + _start != _end)
    {
      assert(false); return false;
    } 
  }
  
  if (psl._qStrand != '-')
  {
    if (psl._qBlockStarts[0] - psl._qChromOffset !=
        _srcStart - psl._qChromOffset)
    {
      assert(false); return false;
    }
    if (psl._qBlockStarts.back() + _blocks.back()._length != psl._qEnd)
    {
      assert(false); return false;
    }
  }
  else
  {
    if (psl._qBlockStarts.back() - psl._qChromOffset !=
        _srcStart - psl._qChromOffset)
    {
      assert(false); return false;
    }
    if (psl._qBlockStarts[0] + _blocks[0]._length != psl._qEnd)
    {
      assert(false); return false;
    }
  }

  return true;
}
