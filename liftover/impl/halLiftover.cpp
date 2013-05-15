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
                       bool addExtraColumns,
                       bool traverseDupes)
{
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _outBedStream = outBedStream;
  _addExtraColumns = addExtraColumns;
  _traverseDupes = traverseDupes;
  _missedSet.clear();
  _tgtSet.clear();
  assert(_srcGenome && inBedStream && tgtGenome && outBedStream);

  _tgtSet.insert(tgtGenome);
  
  if (_inBedVersion == -1)
  {
    _inBedVersion = BedScanner::getBedVersion(inBedStream);
  }
  if (_outBedVersion == -1)
  {
    _outBedVersion = _inBedVersion;
  }

  scan(inBedStream, _inBedVersion);
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
  }
      
  else if (_bedLine._end > (hal_index_t)_srcSequence->getSequenceLength())
  {
    std::cerr << "Skipping interval with endpoint " << _bedLine._end 
              << "because sequence " << _bedLine._chrName << " has length " 
              << _srcSequence->getSequenceLength() << endl;
  }

  liftInterval();
  writeLineResults();
}

void Liftover::writeLineResults()
{

  if (_outBedVersion > 9)
  {
    std::sort(_outBedLines.begin(), _outBedLines.end());
    collapseExtendedBedLines();
  }
  cout << "line " << _lineNumber << endl;  
  for (size_t i = 0; i < _outBedLines.size(); ++i)
  {
    if (_addExtraColumns == false)
    {
      _outBedLines[0]._extra.clear();
    }
    _outBedLines[0].write(*_outBedStream, _outBedVersion);
    _outBedLines[0].write(cout, _outBedVersion);
  }
}

void Liftover::collapseExtendedBedLines()
{

}
