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
               const Genome* refGenome, const Sequence* refSequence,
               hal_index_t refStart, hal_size_t refLength,
               std::set<const Genome*>& targetSet,
               MafExport& mafExport) :
  BedScanner(),
  _mafStream(mafStream),
  _alignment(alignment),
  _refGenome(refGenome),
  _refSequence(refSequence),
  _refStart(refStart),
  _refLength(refLength),
  _targetSet(targetSet),
  _mafExport(mafExport)
{
  if (_refLength == 0)
  {
    if (_refGenome != NULL)
    {
      _refLength = _refGenome->getSequenceLength();
    }
    if (_refSequence != NULL)
    {
      _refLength = _refSequence->getSequenceLength();
    }
  }
}
   
MafBed::~MafBed()
{

}

void MafBed::visitLine()
{
  const Sequence* refSequence = _refGenome->getSequence(_bedLine._chrName);
  if (refSequence != NULL && 
      (_refSequence == NULL || refSequence == _refSequence))
  {
    hal_index_t refStart = _refStart;
    if (_refSequence == NULL)
    {
      // _refStart is in genome coordinate, switch it to be relative to seq
      refStart = _refStart - refSequence->getStartPosition();
    }
    hal_index_t refEnd = std::min(refStart + _refLength, 
                                  refSequence->getSequenceLength());
    if (refStart < 0 || refEnd <= refStart)
    {
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
        hal_index_t start = std::max(_bedLine._start, refStart);
        hal_index_t end = std::min(_bedLine._end, refEnd);
        if (end > start)
        {
          _mafExport.convertSegmentedSequence(_mafStream, _alignment, 
                                              refSequence, start, end - start,
                                              _targetSet);
        }
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
          hal_index_t start = std::max(_bedLine._start +
                                       _bedLine._blocks[i]._start, refStart);
          hal_index_t end = std::min(_bedLine._start +
                                     _bedLine._blocks[i]._start + 
                                     _bedLine._blocks[i]._length, refEnd);
          if (end > start)
          {
            _mafExport.convertSegmentedSequence(_mafStream, _alignment, 
                                                refSequence, start, end - start,
                                                _targetSet);
          }
        }
      }
    }
  }
  else if (_refSequence == NULL)
  {
    cerr << "Line " << _lineNumber << ": BED sequence " << _bedLine._chrName
         << " not found in genome " << _refGenome->getName() << '\n';
  }
}
