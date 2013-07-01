/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "hal.h"
#include "hal4dExtract.h"

using namespace std;
using namespace hal;

Extract4d::Extract4d() 
{

}

Extract4d::~Extract4d()
{

}

void Extract4d::run(const Genome* refGenome, 
                    istream* inBedStream, ostream* outBedStream,
                    int bedVersion)
{
  _refGenome = refGenome;
  _outBedStream = outBedStream;
  _bedVersion = bedVersion;
  if (_bedVersion == -1)
  {
    _bedVersion = getBedVersion(inBedStream);
  }
  scan(inBedStream, _bedVersion);
}

void Extract4d::visitLine()
{
  _outBedLines.clear();
  _refSequence = _refGenome->getSequence(_bedLine._chrName);
  if (_refSequence != NULL)
  {
    if (_bedVersion <= 9)
    {
      if (_bedLine._end <= _bedLine._start ||
          _bedLine._end >= (hal_index_t)_refSequence->getSequenceLength())
      {
        cerr << "Line " << _lineNumber << ": BED coordinates invalid\n";
      }
      else
      {
        extractBed4d();
      }
    }
    else
    {
      extractBlocks4d();
    }
  }
  else if (_refSequence == NULL)
  {
    cerr << "Line " << _lineNumber << ": BED sequence " << _bedLine._chrName
         << " not found in genome " << _refGenome->getName() << '\n';
  }

  write();
}

void Extract4d::extractBed4d()
{
  char c1;
  char c2;
  hal_index_t start = _bedLine._start;
  hal_index_t end = _bedLine._end;
  hal_size_t len = end - start;
  hal_size_t N = len / 3;
  --end;
  if (_bedLine._strand == '-')
  {
    swap(start, end);
  }
  DNAIteratorConstPtr dna = _refSequence->getDNAIterator(start);
  if (_bedLine._strand == '-')
  {
    dna->toReverse();
  }
  
  for (hal_size_t i = 0; i < N; ++i)
  {
    c1 = dna->getChar();
    dna->toRight();
    c2 = dna->getChar();
    dna->toRight();
    if (isFourfoldDegenerate(c1, c2) == true)
    {
      hal_index_t pos = dna->getArrayIndex() - 
         _refSequence->getStartPosition();
      if (_bedLine._strand == '-')
      {
        _outBedLines.push_front(_bedLine);
        _outBedLines.front()._start = pos;
        _outBedLines.front()._end = pos + 1;
      }
      else
      {
        _outBedLines.push_back(_bedLine);
        _outBedLines.back()._start = pos;
        _outBedLines.back()._end = pos + 1;
      }
    }
    dna->toRight();
  }
}

void Extract4d::extractBlocks4d()
{
  _outBedLines.push_back(_bedLine);
  assert(_outBedLines.size() == 1);
  _outBedLines.back()._blocks.clear();

  for (size_t i = 0; i < _bedLine._blocks.size(); ++i)
  {
    if (_bedLine._blocks[i]._length == 0 ||
        _bedLine._blocks[i]._start + _bedLine._blocks[i]._length >= 
        (hal_index_t)_refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
    }
    else
    {
      char c1;
      char c2;
      hal_index_t start = _bedLine._blocks[i]._start;
      hal_index_t end =
         _bedLine._blocks[i]._start + _bedLine._blocks[i]._length;
      hal_size_t N = _bedLine._blocks[i]._length / 3;
      --end;
      if (_bedLine._strand == '-')
      {
        swap(start, end);
      }
      DNAIteratorConstPtr dna = _refSequence->getDNAIterator(start);
      if (_bedLine._strand == '-')
      {
        dna->toReverse();
      }
  
      deque<BedBlock> buffer;
      for (hal_size_t i = 0; i < N; ++i)
      {
        c1 = dna->getChar();
        dna->toRight();
        c2 = dna->getChar();
        dna->toRight();
        if (isFourfoldDegenerate(c1, c2) == true)
        {
          BedBlock block;
          block._start = dna->getArrayIndex() - 
             _refSequence->getStartPosition();
          block._length = 1;
          if (_bedLine._strand == '-')
          {
            buffer.push_front(block);
          }
          else
          {
            buffer.push_back(block);
          }
        }
        dna->toRight();
      }
      for (size_t j = 0; j < buffer.size(); ++j)
      {
        _outBedLines.back()._blocks.push_back(buffer[j]);
      }
    }
  }

  if (_outBedLines.back()._blocks.size() > 0)
  {
    BedLine& bl = _outBedLines.back();
    hal_index_t startDelta = _outBedLines.back()._blocks[0]._start;
    if (startDelta > 0)
    {
      vector<BedBlock>::iterator j;
      vector<BedBlock>::iterator k;
      for (j = bl._blocks.begin(); j != bl._blocks.end(); ++j)
      {
        j->_start -= startDelta;
        k = j;
        ++k;
        assert(k == bl._blocks.end() ||
               k->_start >= (j->_start + j->_length));
      }
      bl._start += startDelta;
    }
    assert(bl._blocks[0]._start == 0);
    
    hal_index_t endDelta = bl._end - (bl._blocks.back()._start + 
                                      bl._blocks.back()._length + 
                                      bl._start);
    assert(endDelta >= 0);
    bl._end -= endDelta;
    assert(bl._blocks.back()._start + 
           bl._blocks.back()._length + bl._start == bl._end);

    assert(_outBedLines.size() == 1);
  }
}

void Extract4d::write()
{
  for (size_t i = 0; i < _outBedLines.size(); ++i)
  {
    _outBedLines[i].write(*_outBedStream, _bedVersion);
  }
}
