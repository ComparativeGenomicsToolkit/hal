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
      throw hal_exception("Only compatible with BED12 input. 4d sites are sensitive to frame."
                          " Even if your genes are all single-exon, please convert them "
                          "to BED12 first.");
    }
    else
    {
      if (_bedLine._end <= _bedLine._start ||
          _bedLine._end > (hal_index_t)_refSequence->getSequenceLength())
      {
        cerr << "Line " << _lineNumber << ": BED coordinates invalid\n";
      }
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

// NB: Throws out 4d sites that occur in a codon split across exon
// boundaries. Otherwise the data would need to be single-copy per
// sequence (not so bad).
void Extract4d::extractConservedBlocks4d()
{
  bool reversed = _bedLine._strand == '-';
  hal_index_t frame = 0;
  deque<BedBlock> buffer;
  bool splitCodon = false;
  // Ugly but way more compact. There is probably a much better way of
  // doing this...
  for (hal_size_t i = reversed ? _bedLine._blocks.size() - 1 : 0;
       reversed ? i != (hal_size_t) -1 : i < _bedLine._blocks.size();
       reversed ? --i : ++i)
  {
    if (_bedLine._blocks[i]._length == 0 ||
        _bedLine._start + _bedLine._blocks[i]._start + _bedLine._blocks[i]._length >
        (hal_index_t)_refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
    }
    else
    {
      hal_index_t start = _bedLine._start + _bedLine._blocks[i]._start;
      hal_index_t end = start + _bedLine._blocks[i]._length;
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
            block._start = colIt->getReferenceSequencePosition() - _bedLine._start;
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

  _outBedLines.push_back(_bedLine);
  assert(_outBedLines.size() == 1);
  _outBedLines.back()._blocks.clear();
  for (size_t j = 0; j < buffer.size(); ++j)
  {
    _outBedLines.back()._blocks.push_back(buffer[j]);
  }

  assert(_outBedLines.size() == 1);
}

// again, keeping this around because of the single-genome alignment
// issue.. otherwise conserving among just 1 genome would be
// equivalent to this
void Extract4d::extractBlocks4d()
{
  hal_index_t frame = 0;
  bool reversed = _bedLine._strand == '-';
  char currCodon[2] = { '\0', '\0' };
  deque<BedBlock> buffer;
  // Ugly but way more compact. There is probably a much better way of
  // doing this...
  for (hal_size_t i = reversed ? _bedLine._blocks.size() - 1 : 0;
       reversed ? i != (hal_size_t) -1 : i < _bedLine._blocks.size();
       reversed ? --i : ++i)
  {
    if (_bedLine._blocks[i]._length == 0 ||
        _bedLine._start + _bedLine._blocks[i]._start + _bedLine._blocks[i]._length >
        (hal_index_t)_refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
    }
    else
    {
      hal_index_t start =  _bedLine._start + _bedLine._blocks[i]._start;
      hal_index_t end = start + _bedLine._blocks[i]._length;
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
              _refSequence->getStartPosition() - _bedLine._start;
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

  _outBedLines.push_back(_bedLine);
  assert(_outBedLines.size() == 1);
  _outBedLines.back()._blocks.clear();
  for (size_t j = 0; j < buffer.size(); ++j)
  {
    _outBedLines.back()._blocks.push_back(buffer[j]);
  }

  assert(_outBedLines.size() == 1);
}

void Extract4d::write()
{
  for (size_t i = 0; i < _outBedLines.size(); ++i)
  {
    _outBedLines[i].write(*_outBedStream, _bedVersion);
  }
}
