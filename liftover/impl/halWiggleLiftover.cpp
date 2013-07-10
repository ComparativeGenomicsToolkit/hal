/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halWiggleLiftover.h"
#include "halBlockMapper.h"

using namespace std;
using namespace hal;

const double WiggleLiftover::DefaultValue = 0.0;
const hal_size_t WiggleLiftover::DefaultTileSize = 10000;

WiggleLiftover::WiggleLiftover()
{

}

WiggleLiftover::~WiggleLiftover()
{

}

void WiggleLiftover::convert(AlignmentConstPtr alignment,
                             const Genome* srcGenome,
                             istream* inputFile,
                             const Genome* tgtGenome,
                             ostream* outputFile,
                             bool traverseDupes,
                             bool unique)
{
  _alignment = alignment;
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _outStream = outputFile;
  _traverseDupes = traverseDupes;
  _unique = unique;
  _srcSequence = NULL;

  if (_srcGenome->getNumTopSegments() > 0)
  {
    _segment = _srcGenome->getTopSegmentIterator();
    _lastIndex = (hal_index_t)_srcGenome->getNumTopSegments();
  }
  else
  {
    _segment = _srcGenome->getBottomSegmentIterator();
    _lastIndex = (hal_index_t)_srcGenome->getNumBottomSegments();
  }

  set<const Genome*> inputSet;
  inputSet.insert(_srcGenome);
  inputSet.insert(_tgtGenome);
  getGenomesInSpanningTree(inputSet, _tgtSet);
  _outVals.init(srcGenome->getSequenceLength(), DefaultValue, DefaultTileSize);
  scan(inputFile);
}

void WiggleLiftover::visitHeader()
{
  _srcSequence = _srcGenome->getSequence(_sequenceName);
  if (_srcSequence == NULL)
  {
    stringstream ss;
    ss << "Sequence " << _sequenceName << " not found in genome " 
       << _srcGenome->getName();
    throw hal_exception(ss.str());
  }
}

void WiggleLiftover::visitLine()
{
  if (_srcSequence == NULL)
  {
    throw hal_exception("Missing Wig header");
  }
  hal_index_t absFirst = _first + _srcSequence->getStartPosition();
  hal_index_t absLast = _last + _srcSequence->getStartPosition();
  if (absFirst < _segment->getStartPosition() || 
      absLast > _segment->getStartPosition())
  {
    mapSegment();
  }
  CoordVal cv = {absFirst, absLast, _value};
  _cvals.push_back(cv);
}
               
void WiggleLiftover::visitEOF()
{
  mapSegment();
}       

void WiggleLiftover::mapSegment()
{
  if (_cvals.empty())
  {
    return;
  }
  
  while (_cvals[0]._first < _segment->getStartPosition())
  {
    assert(_segment->getArrayIndex() > 0);
    assert(_segment->getReversed() == false);
    _segment->toLeft(_cvals[0]._first);
  }

  while (_segment->getArrayIndex() < _lastIndex &&
         _segment->getStartPosition() <= (_cvals.back()._last))
  {
    _segment->getMappedSegments(_mappedSegments, _tgtGenome, &_tgtSet,
                                _traverseDupes);  
    _segment->toRight(_cvals.back()._last);
  }

  vector<MappedSegmentConstPtr> fragments;
  BlockMapper::MSSet emptySet;
  set<hal_index_t> queryCutSet;
  set<hal_index_t> targetCutSet;
  
  for (std::set<MappedSegmentConstPtr>::iterator i = _mappedSegments.begin();
       i != _mappedSegments.end(); ++i)
  {
    BlockMapper::extractSegment(i, emptySet, fragments, &_mappedSegments, 
                                targetCutSet, queryCutSet);
    
  }

}
