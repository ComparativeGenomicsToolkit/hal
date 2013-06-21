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
                       bool traverseDupes)
{
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _outBedStream = outBedStream;
  _addExtraColumns = addExtraColumns;
  _inBedVersion = inBedVersion;
  _outBedVersion = outBedVersion;
  _traverseDupes = traverseDupes;
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
    skipWhiteSpaces(inBedStream);
    std::getline(*inBedStream, firstLineBuffer);
    firstLineStream = new stringstream(firstLineBuffer);
    _inBedVersion = BedScanner::getBedVersion(firstLineStream);
    assert(inBedStream->eof() || _inBedVersion >= 3);
    size_t numCols = BedScanner::getNumColumns(firstLineBuffer);
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
  if (_inBedVersion <= 9 && _outBedVersion > _inBedVersion)
  {
    stringstream ss;
    ss << "Unable to convert from BED version " << _inBedVersion << " to "
       << _outBedVersion << ": input version must be at least as high as output"
       " version.";
    throw hal_exception(ss.str());
  }

  if (firstLineStream != NULL)
  {
    scan(firstLineStream, _inBedVersion);
    delete firstLineStream;
  }
  scan(inBedStream, _inBedVersion);
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

  liftInterval(_outBedLines);  
  
  if (_inBedVersion > 9 && !_bedLine._blocks.empty())
  {
    _mappedBlocks.clear();
    liftBlockIntervals();
    if (_mappedBlocks.size() > 0 && _outBedVersion > 9)
    {
      // extend the intervals
      mergeIntervals();
      // fill them with mapped blocks
      assignBlocksToIntervals();
    }
    if (_outBedVersion <= 9)
    {
      //only map the blocks and forget about the intervals
      writeBlocksAsIntervals();
    }
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
    i->write(*_outBedStream, _outBedVersion);
  }
}

void Liftover::assignBlocksToIntervals()
{
  // sort the mapped blocks by source coordinate
  _mappedBlocks.sort(BedLineSrcLess());

  // sort intervals in target coordinates
  set<BedLine*, BedLinePLess> intervalSet;
  set<BedLine*, BedLinePLess>::iterator setIt;
  for (BedList::iterator i = _outBedLines.begin(); i != _outBedLines.end(); ++i)
  {
    assert((*i)._blocks.empty());
    intervalSet.insert(&*i);
  }
  
  BedList::iterator blockPrev = _mappedBlocks.end();
  BedList::iterator blockNext;
  for (BedList::iterator blockIt = _mappedBlocks.begin(); 
       blockIt != _mappedBlocks.end(); blockIt = blockNext)
  {
    blockNext = blockIt;
    ++blockNext;
    
    // find first interval with start coordinate less than 
    setIt = intervalSet.lower_bound(&*blockIt);    
    if (setIt != intervalSet.begin() && ( 
          setIt == intervalSet.end() ||
          (*setIt)->_chrName != (*blockIt)._chrName ||
          (*setIt)->_start > (*blockIt)._start))
    {
      --setIt;
    }
    assert((*setIt)->_chrName == (*blockIt)._chrName);

    // check if the found interval contains the block
    if ((*setIt)->_start <= (*blockIt)._start &&
        (*setIt)->_end >= (*blockIt)._end)
    {
      bool compatible = blockPrev == _mappedBlocks.end() || 
         (*setIt)->_blocks.empty();         
      if (!compatible)
      {
        hal_index_t blockStart = (*setIt)->_blocks.back()._start + 
           (*setIt)->_start;
        hal_index_t blockEnd = blockStart + (*setIt)->_blocks.back()._length;
        compatible = ((*blockPrev)._start == blockStart && 
                      (*blockPrev)._end == blockEnd &&
                      (*blockIt)._start >= blockEnd) ;
      }
      BedBlock block;
      block._start = (*blockIt)._start - (*setIt)->_start;
      block._length = (*blockIt)._end - (*blockIt)._start;
      
      // we can add the block to the intervals list without breaking
      // ordering in the original input
      if (compatible == true)
      {
        (*setIt)->_blocks.push_back(block);
      }
      else
      {
/*
  cout << "setIt.start " << (*setIt)->_start 
  << " setIt.end " << (*setIt)->_end
  << " blockit.start " << (*blockIt)._start
  << " blockit.end " << (*blockIt)._end << endl;
  cout << " blockprev.start " << (*blockPrev)._start
  << " blockprev.end " << (*blockPrev)._end << endl;
  cout << " blockback " <<( (*setIt)->_blocks.back()._start + 
  (*setIt)->_start)
  << " blockbakcend " 
  << ( (*setIt)->_blocks.back()._start + 
  (*setIt)->_start + (*setIt)->_blocks.back()._length) 
  << endl;
*/
        assert((*setIt)->_blocks.size() > 0);
        // otherwise, we duplicate the containing interval, zap all its
        // blocks, and add the new block.  the old interval can no longer
        // be modified and it is removed from the set
        BedLine newInterval = **setIt;
        intervalSet.erase(setIt);
        newInterval._blocks.clear();
        newInterval._blocks.push_back(block);
        _outBedLines.push_back(newInterval);
        intervalSet.insert(&_outBedLines.back());
      }
    }
    else
    {
      /* cout << "setIt.start " << (*setIt)->_start 
         << " setIt.end " << (*setIt)->_end
         << " blockit.start " << (*blockIt)._start
         << " blockit.end " << (*blockIt)._end << endl;
      */
      assert(false);
    }
    blockPrev = blockIt;
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

// merge intervals (that will contain blocks) if they are on same strand
// and sequence and are colinear (but not necessarily directly adjacent)
void Liftover::mergeIntervals()
{
  assert(_outBedLines.empty() == false);
  if (_outBedLines.size() > 1)
  {
    // sort by target coordinate
    _outBedLines.sort(BedLineLess());

    BedList::iterator i;
    BedList::iterator j;
    for (i = _outBedLines.begin(); i != _outBedLines.end(); ++i)
    {
      j = i;
      ++j;
      if (j != _outBedLines.end())
      {
        if ((*i)._chrName == (*j)._chrName &&
            (*i)._strand == (*j)._strand)
        {
          assert((*i)._start < (*j)._start);
          (*i)._end = max((*i)._end, (*j)._end);
          _outBedLines.erase(j);
        }
      }
    }
  }
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
        if (_bedLine._thickStart != _bedLine._start ||
            _bedLine._thickEnd != _bedLine._end)
        {
          cerr << "Input BED line " << _lineNumber << " warning: "
               << "thickStart different from chromStart or thickEnd "
               << "different from chromEnd not supported.  Assuming they"
               << " are the same." << endl;
        }
      }
      else
      {
        assert(i->_thickStart == 0 && i->_thickEnd == 0);
      }
      
      if (_outBedVersion > 9)
      {
        if  (i->_blocks.size() > 0)
        {
          hal_index_t startDelta = i->_blocks[0]._start;
          assert(startDelta >= 0);
          if (startDelta > 0)
          {
            vector<BedBlock>::iterator j;
            vector<BedBlock>::iterator k;
            for (j = i->_blocks.begin(); j != i->_blocks.end(); ++j)
            {
              j->_start -= startDelta;
              k = j;
              ++k;
              assert(k == i->_blocks.end() ||
                     k->_start >= (j->_start + j->_length));
            }
            i->_start += startDelta;
          }
          assert(i->_blocks[0]._start == 0);

          hal_index_t endDelta = i->_end - (i->_blocks.back()._start + 
                                            i->_blocks.back()._length + 
                                            i->_start);
          assert(endDelta >= 0);
          i->_end -= endDelta;
          assert(i->_blocks.back()._start + 
                 i->_blocks.back()._length + i->_start == i->_end);
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
