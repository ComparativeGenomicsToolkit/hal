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
                    int bedVersion, bool conserved)
{
  _refGenome = refGenome;
  _outBedStream = outBedStream;
  _bedVersion = bedVersion;
  if (_bedVersion == -1)
  {
    _bedVersion = getBedVersion(inBedStream);
  }
  _conserved = conserved;
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
          _bedLine._end > (hal_index_t)_refSequence->getSequenceLength())
      {
        cerr << "Line " << _lineNumber << ": BED coordinates invalid\n";
      }
      else
      {
        if (_conserved) {
          extractConservedBed4d();
        } else {
          extractBed4d();
        }
      }
    }
    else
    {
      if (_conserved) {
        extractConservedBlocks4d();
      } else {
        extractBlocks4d();
      }
    }
  }
  else if (_refSequence == NULL)
  {
    cerr << "Line " << _lineNumber << ": BED sequence " << _bedLine._chrName
         << " not found in genome " << _refGenome->getName() << '\n';
  }

  write();
}

// why keep a unconserved version around?: the column iterator seems
// to fail in the case of a single root genome. Not a common use case
// but it is used in the tests.
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

void Extract4d::extractConservedBed4d()
{
  hal_index_t start = _bedLine._start;
  hal_index_t end = _bedLine._end;
  hal_size_t len = end - start;
  hal_size_t N = len/3;
  --end;
  ColumnIteratorConstPtr colIt = _refSequence->getColumnIterator(NULL, 0, start, NULL_INDEX, false, true, _bedLine._strand == '-');
  
  for (hal_size_t i = 0; i < N; ++i)
  {
    bool is4dSite = true;
    if (_bedLine._strand == '+') {
      // shift to 4d site location
      colIt->toRight();
      colIt->toRight();
    }
    const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
    for (ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
         colMapIt != colMap->end(); ++colMapIt) {
      const ColumnIterator::DNASet *dnaSet = colMapIt->second;
      for (hal_size_t j = 0; j < dnaSet->size(); j++) {
        char c1, c2;
        DNAIteratorConstPtr dna = dnaSet->at(j);
        if ((dna->getReversed() && dna->getArrayIndex() > colMapIt->first->getEndPosition() - 2) ||
            (!dna->getReversed() && dna->getArrayIndex() < colMapIt->first->getStartPosition() + 2)) {
          is4dSite = false;
          break;
        }
        dna->toLeft();
        c2 = dna->getChar();
        dna->toLeft();
        c1 = dna->getChar();
        if (!isFourfoldDegenerate(c1, c2)) {
          is4dSite = false;
          break;
        }
      }
    }
    if (is4dSite)
    {
      hal_index_t pos = colIt->getReferenceSequencePosition();
      _outBedLines.push_back(_bedLine);
      _outBedLines.back()._start = pos;
      _outBedLines.back()._end = pos + 1;
    }
    if (_bedLine._strand == '-') {
      colIt->toRight();
      colIt->toRight();
    }
    if (!colIt->lastColumn()) {
      colIt->toRight();
    }
  }
}

// NB: Throws out 4d sites that occur in a codon split across exon
// boundaries. Otherwise the data would need to be single-copy per
// sequence (not so bad).
void Extract4d::extractConservedBlocks4d()
{
  _outBedLines.push_back(_bedLine);
  assert(_outBedLines.size() == 1);
  _outBedLines.back()._blocks.clear();
  bool reversed = _bedLine._strand == '-';
  hal_index_t frame = 0;
  deque<BedBlock> buffer;
  bool splitCodon = false;
  // Ugly but way more compact. There is probably a much better way of
  // doing this...
  for (int64_t i = reversed ? _bedLine._blocks.size() - 1 : 0;
       reversed ? i >= 0 : i < _bedLine._blocks.size();
       reversed ? --i : ++i)
  {
    if (_bedLine._blocks[i]._length == 0 ||
        _bedLine._blocks[i]._start + _bedLine._blocks[i]._length > 
        (hal_index_t)_refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
    }
    else
    {
      hal_index_t start = _bedLine._blocks[i]._start;
      hal_index_t end =
        _bedLine._blocks[i]._start + _bedLine._blocks[i]._length;
      --end;
      if (reversed) {
        swap(start, end);
      }
      ColumnIteratorConstPtr colIt = _refSequence->getColumnIterator(NULL, 0, start, NULL_INDEX, false, true, _bedLine._strand == '-');
      for (hal_index_t n = 0; n < _bedLine._blocks[i]._length; ++n)
      {
        if (frame == 2) {
          bool is4dSite = true;
          const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
          for (ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
               colMapIt != colMap->end(); ++colMapIt) {
            const ColumnIterator::DNASet *dnaSet = colMapIt->second;
            for (hal_size_t j = 0; j < dnaSet->size(); j++) {
              char c1, c2;
              DNAIteratorConstPtr dna = dnaSet->at(j);
              if ((dna->getReversed() && dna->getArrayIndex() > colMapIt->first->getEndPosition() - 2) ||
                  (!dna->getReversed() && dna->getArrayIndex() < colMapIt->first->getStartPosition() + 2)) {
                is4dSite = false;
                break;
              }
              dna->toLeft();
              c2 = dna->getChar();
              dna->toLeft();
              c1 = dna->getChar();
              if (!isFourfoldDegenerate(c1, c2)) {
                is4dSite = false;
                break;
              }
            }
            frame = 0;
          }
          if (is4dSite && !splitCodon)
          {
            BedBlock block;
            block._start = colIt->getReferenceSequencePosition();
            block._length = 1;
            if (reversed) {
              buffer.push_front(block);
            } else {
              buffer.push_back(block);
            }
          }
        } else {
          frame++;
        }
        splitCodon = false;
        if (n != _bedLine._blocks[i]._length - 1) {
          if (reversed) {
            colIt->toSite(colIt->getReferenceSequencePosition() + _refSequence->getStartPosition() - 1,
                          colIt->getReferenceSequencePosition() + _refSequence->getStartPosition(),
                          true);
          } else {
            colIt->toSite(colIt->getReferenceSequencePosition() + _refSequence->getStartPosition() + 1,
                          colIt->getReferenceSequencePosition() + _refSequence->getStartPosition() + 2,
                          true);
          }
        }
      }
      if (frame != 0) {
        splitCodon = true;
      }
    }
  }
  for (size_t j = 0; j < buffer.size(); ++j)
  {
    _outBedLines.back()._blocks.push_back(buffer[j]);
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

// again, keeping this around because of the single-genome alignment
// issue.. otherwise conserving among just 1 genome would be
// equivalent to this
void Extract4d::extractBlocks4d()
{
  _outBedLines.push_back(_bedLine);
  assert(_outBedLines.size() == 1);
  _outBedLines.back()._blocks.clear();
  hal_index_t frame = 0;
  bool reversed = _bedLine._strand == '-';
  char currCodon[2] = { '\0', '\0' };
  deque<BedBlock> buffer;
  // Ugly but way more compact. There is probably a much better way of
  // doing this...
  for (int64_t i = reversed ? _bedLine._blocks.size() - 1 : 0;
       reversed ? i >= 0 : i < _bedLine._blocks.size();
       reversed ? --i : ++i)
  {
    if (_bedLine._blocks[i]._length == 0 ||
        _bedLine._blocks[i]._start + _bedLine._blocks[i]._length > 
        (hal_index_t)_refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
    }
    else
    {
      hal_index_t start = _bedLine._blocks[i]._start;
      hal_index_t end =
         _bedLine._blocks[i]._start + _bedLine._blocks[i]._length;
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
  
      for (hal_index_t n = 0; n < _bedLine._blocks[i]._length; ++n)
      {
        if (frame == 2) {
          if (isFourfoldDegenerate(currCodon[0], currCodon[1]) == true)
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
          frame = 0;
        } else {
          currCodon[frame] = dna->getChar();
          frame++;
        }
        dna->toRight();
      }
    }
  }

  for (size_t j = 0; j < buffer.size(); ++j)
  {
    _outBedLines.back()._blocks.push_back(buffer[j]);
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
