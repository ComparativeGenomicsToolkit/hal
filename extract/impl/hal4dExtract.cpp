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
        extractBlocks4d(true);
      } else {
        extractBlocks4d(false);
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

// Check if the 4d site is still a 4d site in all aligned regions.
static bool is4dSiteConserved(const Sequence *sequence, hal_index_t pos, bool reversed)
{
  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(NULL, 0, pos, NULL_INDEX, false, true, reversed);
  bool isConserved = true;
  const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
  for (ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
       colMapIt != colMap->end(); ++colMapIt) {
    const ColumnIterator::DNASet *dnaSet = colMapIt->second;
    for (hal_size_t j = 0; j < dnaSet->size(); j++) {
      char c1, c2;
      DNAIteratorConstPtr dna = dnaSet->at(j);
      if ((dna->getReversed() && dna->getArrayIndex() > colMapIt->first->getEndPosition() - 2) ||
          (!dna->getReversed() && dna->getArrayIndex() < colMapIt->first->getStartPosition() + 2)) {
        isConserved = false;
        break;
      }
      dna->toLeft();
      c2 = dna->getChar();
      dna->toLeft();
      c1 = dna->getChar();
      if (!isFourfoldDegenerate(c1, c2)) {
        isConserved = false;
        break;
      }
    }
  }
  return isConserved;
}

// NB: If conserved == true, throws out 4d sites that occur in a codon
// split across exon boundaries. Otherwise the data would need to be
// single-copy per sequence (not so bad).
void Extract4d::extractBlocks4d(bool conserved)
{
  hal_index_t frame = 0;
  hal_index_t cdsStart = _bedLine._thickStart;
  hal_index_t cdsEnd = _bedLine._thickEnd;
  bool reversed = _bedLine._strand == '-';
  char currCodonPrefix[2] = { '\0', '\0' };
  deque<BedBlock> buffer;
  for (hal_size_t i = reversed ? _bedLine._blocks.size() - 1 : 0;
       reversed ? i != (hal_size_t) -1 : i < _bedLine._blocks.size();
       reversed ? --i : ++i)
  {
    BedBlock block = _bedLine._blocks[i];
    if (block._length == 0 ||
        _bedLine._start + block._start + block._length >
        (hal_index_t)_refSequence->getSequenceLength())
    {
      // Bad bed block.
      stringstream ss;
      ss << "Line " << _lineNumber << ", block " << i 
           <<": BED coordinates invalid\n";
      throw hal_exception(ss.str());
    }
    hal_index_t start =  _bedLine._start + block._start;
    hal_index_t end = start + block._length;
    // Force start and end to be within the CDS start/end.
    if (end <= cdsStart) {
      continue;
    } else {
      start = max(cdsStart, start);
    }
    if (start >= cdsEnd) {
      continue;
    } else {
      end = min(end, cdsEnd);
    }
    hal_index_t length = end - start;
    if (reversed)
    {
      --end;
      swap(start, end);
    }
    DNAIteratorConstPtr dna = _refSequence->getDNAIterator(start);
    if (reversed)
    {
      dna->toReverse();
    }

    for (hal_index_t n = 0; n < length; ++n)
    {
      if (frame == 2) {
        if (isFourfoldDegenerate(currCodonPrefix[0], currCodonPrefix[1]))
        {
          if (!conserved ||
              // We can't deal with split codons in conserved mode currently.
              (n >= 2 &&
               // Check if the 4d site is a 4d site in all species.
               is4dSiteConserved(_refSequence, dna->getArrayIndex() - _refSequence->getStartPosition(), reversed)))
          {
            BedBlock outBlock;
            outBlock._start = dna->getArrayIndex() - 
              _refSequence->getStartPosition() - _bedLine._start;
            outBlock._length = 1;
            if (reversed)
            {
              buffer.push_front(outBlock);
            }
            else
            {
              buffer.push_back(outBlock);
            }
          }
        }
        frame = 0;
      } else {
        currCodonPrefix[frame] = dna->getChar();
        frame++;
      }
      dna->toRight();
    }
  }

  if (buffer.size() > 0) {
    // We have at least one block to write, put them into an output bed line.
    _outBedLines.push_back(_bedLine);
    assert(_outBedLines.size() == 1);
    _outBedLines.back()._blocks.clear();
    for (size_t j = 0; j < buffer.size(); ++j)
    {
      _outBedLines.back()._blocks.push_back(buffer[j]);
    }

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
