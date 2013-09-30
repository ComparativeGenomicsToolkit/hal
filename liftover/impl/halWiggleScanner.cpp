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
#include <cctype>

#include "halWiggleScanner.h"

using namespace std;
using namespace hal;

WiggleScanner::WiggleScanner() : _wiggleStream(NULL)
{

}

WiggleScanner::~WiggleScanner()
{

}

void WiggleScanner::scan(const string& wigglePath)
{
  assert(_wiggleStream == NULL);
  _wiggleStream = new ifstream(wigglePath.c_str());
  
  try {
    scan(_wiggleStream);
  }
  catch(hal_exception e)
  {
    delete _wiggleStream;
    _wiggleStream = NULL;
    stringstream ss;
    ss << e.what() << " in file " << wigglePath;
    throw hal_exception(ss.str());
  }  

  delete _wiggleStream;
  _wiggleStream = NULL;
}

void WiggleScanner::scan(istream* is)
{
  visitBegin();
  _wiggleStream = is;

  if (_wiggleStream->bad())
  {
    throw hal_exception("Error reading wiggle input stream");
  }
  string lineBuffer;
  _lineNumber = 0;
  try
  {
    skipWhiteSpaces(_wiggleStream);
    while (_wiggleStream->good())
    {
      ++_lineNumber;
      std::getline(*is, lineBuffer);
      if (scanHeader(lineBuffer) == true)
      {
        visitHeader();
      }
      else
      {
        scanLine(lineBuffer);
        visitLine();
      }
      skipWhiteSpaces(_wiggleStream);
    }
  }
  catch(hal_exception e)
  {
    stringstream ss;
    ss << e.what() << " -- input wiggle line " << _lineNumber;
    throw hal_exception(ss.str());
  }
  visitEOF();
  _wiggleStream = NULL;
}

bool WiggleScanner::scanHeader(const string& lineBuffer)
{
  stringstream ss(lineBuffer);
  ss >> _buffer;
  if (ss.good() && _buffer == "variableStep")
  {
    _fixedStep = false;
    ss >> _buffer;
    if (!ss || _buffer.length() <= 6 || _buffer.substr(0, 6) != "chrom=")
    {
      throw hal_exception("Error parsing chrom in variableStep header");
    }
    _sequenceName = _buffer.substr(6);
    
    ss >> _buffer;
    if (!ss || _buffer.length() <= 5 || _buffer.substr(0, 5) != "span=")
    {
      _span = NULL_INDEX;
    }
    else
    {
      _buffer = _buffer.substr(5);
      stringstream ss1(_buffer);
      ss1 >> _span;
      if (!ss1)
      {
        _span = NULL_INDEX;
      }
    }
    return true;
  }
  else if (ss.good() && _buffer == "fixedStep")
  {
    _fixedStep = true;
    _offset = 0;
    ss >> _buffer;
    if (!ss || _buffer.length() <= 6 || _buffer.substr(0, 6) != "chrom=")
    {
      throw hal_exception("Error parsing chrom in fixedStep header");
    }
    _sequenceName = _buffer.substr(6);
    
    ss >> _buffer;
    if (!ss || _buffer.length() <= 6 || _buffer.substr(0, 6) != "start=")
    {
      throw hal_exception("Error parsing start in fixedStep header");
    }
    _buffer = _buffer.substr(6);
    stringstream ss1(_buffer);
    ss1 >> _start;
    if (!ss1)
    {
      throw hal_exception("Error parsing start in fixedStep header");
    }
    assert(_start > 0);
    // store internally in 0-based coordinates.
    --_start;

    ss >> _buffer;
    if (!ss || _buffer.length() <= 5 || _buffer.substr(0, 5) != "step=")
    {
      throw hal_exception("Error parsing step in fixedStep header");
    }
    _buffer = _buffer.substr(5);
    stringstream ss2(_buffer);
    ss2 >> _step;
    if (!ss2)
    {
      throw hal_exception("Error parsing step in fixedStep header");
    }

    ss >> _buffer;
    if (!ss || _buffer.length() <= 5 || _buffer.substr(0, 5) != "span=")
    {
      _span = NULL_INDEX;
    }
    else
    {
      _buffer = _buffer.substr(5);
      stringstream ss1(_buffer);
      ss1 >> _span;
      if (!ss1)
      {
        _span = NULL_INDEX;
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

void WiggleScanner::scanLine(const string& lineBuffer)
{
  stringstream ss(lineBuffer);
  if (_fixedStep == true)
  {
    _first = _start + _offset * _step;
    ++_offset;
  }
  else
  {
    ss >> _first;
    if (!ss)
    {
      stringstream ss1;
      ss1 << "Error parsing position for " << _sequenceName;
      throw hal_exception(ss1.str());
    }
    assert(_first > 0);
    // store internally in 0-based coordinates.
    --_start;
  }

  ss >> _value;
  if (!ss)
  {
    stringstream ss1;
    ss1 << "Error parsing value for " << _sequenceName << " pos " << _start;
    throw hal_exception(ss1.str());
  }
  _last = _first;
  if (_span > 1)
  {
    _last += _span - 1;
  }
}

void WiggleScanner::visitBegin()
{
}

void WiggleScanner::visitLine()
{
}

void WiggleScanner::visitHeader()
{
}

void WiggleScanner::visitEOF()
{
}

void WiggleScanner::skipWhiteSpaces(istream* wiggleStream)
{
  while (wiggleStream->good() && std::isspace(wiggleStream->peek()))
  {
    wiggleStream->get();
  }
}
