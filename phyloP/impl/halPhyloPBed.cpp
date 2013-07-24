/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <algorithm>
#include "halPhyloPBed.h"
// get rid of phast macros in order to use stl
#undef min
#undef max
using namespace std;
using namespace hal;


PhyloPBed::PhyloPBed(AlignmentConstPtr alignment,
                     const Genome* refGenome, const Sequence* refSequence,
                     hal_index_t start, hal_size_t length, 
                     hal_size_t step, PhyloP& phyloP, ostream& outStream) :
  BedScanner(),
  _alignment(alignment),
  _refGenome(refGenome),
  _refSequence(refSequence),
  _refStart(start),
  _refLength(length),
  _step(step),
  _phyloP(phyloP),
  _outStream(outStream)
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
   
PhyloPBed::~PhyloPBed()
{

}

void PhyloPBed::visitLine()
{
  const Sequence* refSequence = _refGenome->getSequence(_bedLine._chrName);
  if (refSequence != NULL && 
      (_refSequence == NULL || refSequence == _refSequence))
  {
    hal_index_t refStart = _refStart;
    if (_refSequence == NULL)
    {
      // _refStart is in genome coordinate, switch it to be relative to sequence
      refStart = _refStart - refSequence->getStartPosition();
    }
    hal_index_t refEnd = std::min(refStart + _refLength, 
                                  (hal_index_t)
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
          _phyloP.processSequence(refSequence, start, end - start, _step);
        }
      }
    }
    else
    {
      for (size_t i = 0; i < _bedLine._blocks.size(); ++i)
      {
        if (_bedLine._blocks[i]._length == 0 ||
            _bedLine._blocks[i]._start + _bedLine._blocks[i]._length >= 
            (hal_index_t)refSequence->getSequenceLength())
        {
          cerr << "Line " << _lineNumber << ", block " << i 
               <<": BED coordinates invalid\n";
        }
        else
        {
          hal_index_t start = std::max(_bedLine._blocks[i]._start, refStart);
          hal_index_t end = std::min(_bedLine._blocks[i]._start + 
                                     _bedLine._blocks[i]._length, refEnd);
          if (end > start)
          {
            _phyloP.processSequence(refSequence, start, end - start, _step);
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
