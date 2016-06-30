/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halMafBed.h"

using namespace std;
using namespace hal;

MafBed::MafBed(std::ostream& mafStream, AlignmentConstPtr alignment,
               const Genome* refGenome, std::set<const Genome*>& targetSet,
               MafExport& mafExport) :
  BedScanner(),
  _mafStream(mafStream),
  _alignment(alignment),
  _refGenome(refGenome),
  _targetSet(targetSet),
  _mafExport(mafExport)
{
}
   
MafBed::~MafBed()
{

}

void MafBed::visitLine()
{
  const Sequence* refSequence = _refGenome->getSequence(_bedLine._chrName);
  if (refSequence == NULL)
  {
    cerr << "Line " << _lineNumber << ": BED sequence " << _bedLine._chrName
         << " not found in genome " << _refGenome->getName() << '\n';
    return;
  }
  if (_bedVersion <= 9)
  {
    if (_bedLine._end <= _bedLine._start ||
        _bedLine._end > (hal_index_t)refSequence->getSequenceLength())
    {
      cerr << "Line " << _lineNumber << ": BED coordinates invalid\n";
    }
    else
    {
      hal_index_t start = _bedLine._start;
      hal_index_t end = _bedLine._end;
      _mafExport.convertSegmentedSequence(_mafStream, _alignment, 
                                          refSequence, start, end - start,
                                          _targetSet);
    }
  }
  else
  {
    for (size_t i = 0; i < _bedLine._blocks.size(); ++i)
    {
      if (_bedLine._blocks[i]._length == 0 ||
          _bedLine._start + _bedLine._blocks[i]._start +
          _bedLine._blocks[i]._length >= 
          (hal_index_t)refSequence->getSequenceLength())
      {
        cerr << "Line " << _lineNumber << ", block " << i 
             <<": BED coordinates invalid\n";
      }
      else
      {
        hal_index_t start = _bedLine._start +_bedLine._blocks[i]._start;
        hal_index_t end = _bedLine._start + _bedLine._blocks[i]._start + _bedLine._blocks[i]._length;
        _mafExport.convertSegmentedSequence(_mafStream, _alignment, 
                                            refSequence, start, end - start,
                                            _targetSet);
      }
    }
  }
}
