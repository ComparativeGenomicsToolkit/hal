/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halLiftover.h"

using namespace std;
using namespace hal;

Liftover::Liftover() : _outBedStream(NULL),                       
                       _inBedVersion(-1), _outBedVersion(-1),
                       _outPSL(false), _outPSLWithName(false),
                       _srcGenome(NULL), _tgtGenome(NULL)
{

}

Liftover::~Liftover()
{

}

void Liftover::convert(AlignmentConstPtr alignment,
                       const Genome* srcGenome,
                       istream* inBedStream,
                       const Genome* tgtGenome,
                       ostream* outBedStream,
                       int inBedVersion,
                       int outBedVersion,
                       bool addExtraColumns,
                       bool traverseDupes,
                       bool outPSL,
                       bool outPSLWithName,
                       const locale* inLocale,
                       const Genome *coalescenceLimit)
{
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _coalescenceLimit = coalescenceLimit;
  _outBedStream = outBedStream;
  _addExtraColumns = addExtraColumns;
  _inBedVersion = inBedVersion;
  _outBedVersion = outBedVersion;
  _traverseDupes = traverseDupes;
  _outPSL = outPSL;
  _outPSLWithName = outPSLWithName;
  _inLocale = inLocale;
  _missedSet.clear();
  _tgtSet.clear();
  assert(_srcGenome && inBedStream && tgtGenome && outBedStream);

  _tgtSet.insert(tgtGenome);
  
  // we copy into a stringstream because we never want to
  // run getBedVersion (which does random access) on cin (which is a
  // possible value for inBedStream)
  string firstLineBuffer;
  stringstream* firstLineStream = NULL;
  if (_inBedVersion <= 0)
  {
    skipWhiteSpaces(inBedStream, _inLocale);
    std::getline(*inBedStream, firstLineBuffer);
    firstLineStream = new stringstream(firstLineBuffer);
    _inBedVersion = BedScanner::getBedVersion(firstLineStream, inLocale);
    assert(inBedStream->eof() || _inBedVersion >= 3);
    size_t numCols = BedScanner::getNumColumns(firstLineBuffer, inLocale);
    if ((int)numCols > _inBedVersion)
    {
      cerr << "Warning: auto-detecting input BED version " << _inBedVersion
           << " even though " << numCols << " columns present" << endl;
    }
  }
  if (_outBedVersion <= 0)
  {
    _outBedVersion = _inBedVersion;
  }

  if (firstLineStream != NULL)
  {
    scan(firstLineStream, _inBedVersion, inLocale);
    delete firstLineStream;
  }
  scan(inBedStream, _inBedVersion, inLocale);
}

void Liftover::visitBegin()
{
}

void Liftover::visitLine()
{
  _outBedLines.clear();
  _srcSequence = _srcGenome->getSequence(_bedLine._chrName);
  if (_srcSequence == NULL)
  {
    pair<set<string>::iterator, bool> result = _missedSet.insert(
      _bedLine._chrName);
    if (result.second == true)
    {
      std::cerr << "Unable to find sequence " << _bedLine._chrName 
                << " in genome " << _srcGenome->getName() << endl;
    }
    return;
  }
      
  else if (_bedLine._end > (hal_index_t)_srcSequence->getSequenceLength())
  {
    std::cerr << "Skipping interval with endpoint " << _bedLine._end 
              << "because sequence " << _bedLine._chrName << " has length " 
              << _srcSequence->getSequenceLength() << endl;
    return;
  }

  else if (_inBedVersion > 9 && _bedLine._blocks.empty())
  {
    std::cerr << "Skipping input line with 0 blocks" << endl;
    return;
  }

  _mappedBlocks.clear();
  if (_inBedVersion <= 9)
  {
    liftInterval(_mappedBlocks);  
  }
  else
  {
    assert(!_bedLine._blocks.empty());
    liftBlockIntervals();
  }

  if (_mappedBlocks.size() > 0 && _outBedVersion > 9)
  {
    // fill them with mapped blocks
    assignBlocksToIntervals();
  }
  if (_outBedVersion <= 9)
  {
    //only map the blocks and forget about the intervals
    writeBlocksAsIntervals();
  }

  cleanResults();
  _outBedLines.sort(BedLineSrcLess());
  writeLineResults();
}

void Liftover::visitEOF()
{
}

void Liftover::writeLineResults()
{
  BedList::iterator i = _outBedLines.begin();
  for (; i != _outBedLines.end(); ++i)
  {
    if (_addExtraColumns == false)
    {
      i->_extra.clear();
    }
    if (_outPSL == false)
    {
      i->write(*_outBedStream, _outBedVersion);
    }
    else
    {
      i->writePSL(*_outBedStream, _outPSLWithName);
    }
  }
}

void Liftover::assignBlocksToIntervals()
{
  assert(_outBedLines.size() == 0);

   // sort the mapped blocks by source coordinate
  _mappedBlocks.sort(BedLineSrcLess());
  hal_index_t prevSrcBlockEnd = NULL_INDEX;
  BedList::iterator blockNext;
  for (BedList::iterator blockIt = _mappedBlocks.begin();
       blockIt != _mappedBlocks.end(); ++blockIt)
  {
    blockNext = blockIt;
    ++blockNext;
    hal_index_t srcBlockEnd = blockIt->_srcStart +
      (blockIt->_end - blockIt->_start);
    bool dupe = (blockIt->_srcStart < prevSrcBlockEnd)  ||
      ((blockNext != _mappedBlocks.end()) &&
       (blockNext->_srcStart < srcBlockEnd));
    if (_outBedLines.empty() || 
        // filter dupes in psl but let them be on single bed line
        (_outPSL && dupe) ||
        !compatible(_outBedLines.back(), *blockIt))
    {
      _outBedLines.push_back(*blockIt);
    }
    prevSrcBlockEnd = blockIt->_srcStart + (blockIt->_end - blockIt->_start);
    BedLine& tgtBed = _outBedLines.back();
    tgtBed._start = min(tgtBed._start, blockIt->_start);
    tgtBed._end = max(tgtBed._end, blockIt->_end);
    // keep start absolute for now
    BedBlock block;
    block._start = blockIt->_start;
    block._length = blockIt->_end - blockIt->_start;
    tgtBed._blocks.push_back(block);

    if (_outPSL == true)
    {
      assert(tgtBed._psl.size() == 1);
      tgtBed._psl[0]._qBlockStarts.push_back(blockIt->_srcStart);
      // note that these get done on assignmnet for first block
      if (tgtBed._blocks.size() > 1)
      {
        tgtBed._psl[0]._matches += blockIt->_psl[0]._matches;
        tgtBed._psl[0]._misMatches += blockIt->_psl[0]._misMatches;
        tgtBed._psl[0]._repMatches += blockIt->_psl[0]._repMatches;
        tgtBed._psl[0]._nCount += blockIt->_psl[0]._nCount;
      }
      assert(tgtBed._blocks.size() == 
             tgtBed._psl[0]._qBlockStarts.size());
    }
  }

  // relativize block starts
  for (BedList::iterator blockIt = _outBedLines.begin();
       blockIt != _outBedLines.end(); ++blockIt)
  {
    for (size_t i = 0; i < blockIt->_blocks.size(); ++i)
    {
      assert(blockIt->_blocks[i]._start >= blockIt->_start);
      blockIt->_blocks[i]._start -= blockIt->_start;
    }
  }

  // Flip the block ordering to make sure it's ascending in the output
  if (!_outBedLines.empty())
  {
    flipBlocks(_outBedLines);
  }

  // Fill in some informtation for the insert PSL fields
  if (_outPSL == true)
  {
    computePSLInserts(_outBedLines);
  }
}

bool Liftover::compatible(const BedLine& tgtBed, const BedLine& newBlock)
{
  if (tgtBed._strand != newBlock._strand)
  {
    return false;
  }
  assert(newBlock._srcStart >= tgtBed._srcStart);
  if (tgtBed._srcStart == newBlock._srcStart)
  {
    return false;
  }

  hal_index_t delta;
  const BedBlock& tgtBlock = tgtBed._blocks.back();

  if (tgtBed._strand != _bedLine._strand)
  {
    delta = tgtBlock._start - newBlock._end;
  }
  else
  {
    delta = newBlock._start - 
       (tgtBlock._start + (hal_index_t)tgtBlock._length);
  }
  if (delta < 0)
  {
    return false;
  }

  if (tgtBed._chrName != newBlock._chrName)
  {
    return false;
  }

  return true;
}

void Liftover::flipBlocks(BedList& bedList)
{
  for (BedList::iterator bedIt = bedList.begin(); bedIt != bedList.end(); 
       ++bedIt)
  {
    if (bedIt->_blocks.size() > 1)
    {
      hal_index_t delta = bedIt->_blocks[1]._start - 
         (bedIt->_blocks[0]._start + (hal_index_t)bedIt->_blocks[0]._length);
      bool mustFlip = false;

      if (_outPSL != true)
      {
        mustFlip = delta < 0;
      }      
      else
      {
        mustFlip = (bedIt->_strand == '-' && delta >= 0 || 
                    bedIt->_strand != '-' && delta < 0);
      }
      
      if (mustFlip == true)
      {         
        std::reverse(bedIt->_blocks.begin(),
                     bedIt->_blocks.end());
        if (_outPSL == true)
        {
          std::reverse(bedIt->_psl[0]._qBlockStarts.begin(),
                       bedIt->_psl[0]._qBlockStarts.end());
        }
      }

#ifndef NDEBUG
      if (_outPSL == true)
      {
        for (size_t i = 1; i < bedIt->_blocks.size(); ++i)
        {
          if (bedIt->_strand == '-')
          {
            assert(bedIt->_blocks[i]._start < bedIt->_blocks[i-1]._start);
          }
          else
          {
            assert(bedIt->_blocks[i]._start > bedIt->_blocks[i-1]._start);
          }
          if (bedIt->_psl[0]._qStrand == '-')
          {
            assert(bedIt->_psl[0]._qBlockStarts[i] < 
                   bedIt->_psl[0]._qBlockStarts[i-1]);
          }
          else
          {
            assert(bedIt->_psl[0]._qBlockStarts[i] >
                   bedIt->_psl[0]._qBlockStarts[i-1]);
          }
        }
      }
#endif      
    }
  }
}

void Liftover::computePSLInserts(BedList& bedList)
{
  for (BedList::iterator bedIt = bedList.begin(); bedIt != bedList.end(); 
       ++bedIt)
  {
    PSLInfo& psl = bedIt->_psl[0];
    psl._qNumInsert = 0;
    psl._qBaseInsert = 0;
    psl._tNumInsert = 0;
    psl._tBaseInsert = 0;  

    assert(bedIt->_blocks.size() == psl._qBlockStarts.size());
    vector<BedBlock>::iterator blockIt = bedIt->_blocks.begin();
    vector<BedBlock>::iterator blockPrev = blockIt;

    vector<hal_index_t>::iterator qStartIt = psl._qBlockStarts.begin();
    vector<hal_index_t>::iterator qStartPrev = qStartIt;

    if (blockIt != bedIt->_blocks.end())
    {
      ++blockIt;
      ++qStartIt;
    }       
    for (; blockIt != bedIt->_blocks.end(); ++blockIt, ++blockPrev, 
            ++qStartIt, ++qStartPrev)
    {
      if (bedIt->_strand == '-')
      {
        swap(blockIt, blockPrev);
      }
      assert(blockIt->_start >= (blockPrev->_start + blockPrev->_length));
      hal_size_t gap = blockIt->_start -
         (blockPrev->_start + blockPrev->_length);
      if (gap > 0)
      {
        ++psl._tNumInsert;
        psl._tBaseInsert += gap;
      }      
      if (bedIt->_strand == '-')
      {
        swap(blockIt, blockPrev);
      }

      if (psl._qStrand == '-')
      {
        swap(qStartIt, qStartPrev);
        swap(blockIt, blockPrev);
      }

      if (*qStartIt >= (*qStartPrev + blockPrev->_length))
      {
        gap = *qStartIt - (*qStartPrev + blockPrev->_length);
      }
      else
      {
        // Duplicated blocks can overlap.
        gap = 0;
      }
      if (gap > 0)
      {
        ++psl._qNumInsert;
        psl._qBaseInsert += gap;
      }      
      if (psl._qStrand == '-')
      {
        swap(qStartIt, qStartPrev);
        swap(blockIt, blockPrev);
      }
    }
  }
}
  

void Liftover::writeBlocksAsIntervals()
{
  _outBedLines = _mappedBlocks;
}

void Liftover::liftBlockIntervals()
{
  BedLine originalBedLine = _bedLine;
  std::sort(_bedLine._blocks.begin(), _bedLine._blocks.end());
  vector<BedBlock>::iterator blockIt = _bedLine._blocks.begin();
  for (; blockIt != _bedLine._blocks.end(); ++blockIt)
  {
    _bedLine._start = blockIt->_start + originalBedLine._start;
    _bedLine._end = _bedLine._start + blockIt->_length;
    if (_bedLine._end > _bedLine._start)
    {
      liftInterval(_mappedBlocks);
    }
  }
  _bedLine._start = originalBedLine._start;
  _bedLine._end = originalBedLine._end;
}

// do a little postprocessing on the lifted intervals to make sure
// they are bed compliant
void Liftover::cleanResults()
{
  if (_outBedVersion > 6)
  {
    BedList::iterator i;
    BedList::iterator j;
    for (i = _outBedLines.begin(); i != _outBedLines.end(); i = j)
    {
      j = i;
      ++j;
      if (_bedLine._thickStart != 0 || _bedLine._thickEnd != 0)
      {
        i->_thickStart = i->_start;
        i->_thickEnd = i->_end;
/*        if ((_bedLine._thickStart != _bedLine._start ||
            _bedLine._thickEnd != _bedLine._end) && !_warnedThickStart)
        {
          cerr << "Input BED line " << _lineNumber << " warning: "
               << "thickStart different from chromStart or thickEnd "
               << "different from chromEnd not supported.  Assuming they"
               << " are the same." << endl;
               }*/
      }
      else
      {
        assert(i->_thickStart == 0 && i->_thickEnd == 0);
      }
      
      if (_outBedVersion > 9)
      {
        if  (i->_blocks.size() > 0)
        {
          if (_outPSL == true)
          {
            assert(i->_psl.size() == 1);
            i->_srcStart = numeric_limits<hal_index_t>::max();
            i->_psl[0]._qEnd = 0;
            for (size_t j = 0; j < i->_psl[0]._qBlockStarts.size(); ++j)
            {
              i->_srcStart = min(i->_srcStart, i->_psl[0]._qBlockStarts[j]);
              i->_psl[0]._qEnd = max(i->_psl[0]._qEnd, 
                                     (hal_size_t)i->_psl[0]._qBlockStarts[j] +
                                     i->_blocks[j]._length);
            }
          }
        }
        else
        {
          // if we're in bed 12, we don't want any empty regions in the 
          // output
          _outBedLines.erase(i);
        }
      }
    }
  }
}
