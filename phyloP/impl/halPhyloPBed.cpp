/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halPhyloPBed.h"

using namespace std;
using namespace hal;


PhyloPBed::PhyloPBed(AlignmentConstPtr alignment,
                     const Genome* refGenome, const Sequence* refSequence,
                     hal_size_t step, PhyloP& phyloP, ostream& outStream) :
  BedScanner(),
  _alignment(alignment),
  _refGenome(refGenome),
  _refSequence(refSequence),
  _step(step),
  _phyloP(phyloP),
  _outStream(outStream)
{
  
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
    if (_bedVersion <= 9)
    {
      if (_bedLine._end <= _bedLine._start ||
          _bedLine._end > (hal_index_t)refSequence->getSequenceLength())
      {
        cerr << "Line " << _lineNumber << ": BED coordinates invalid\n";
      }
      else
      {
        _phyloP.processSequence(refSequence, _bedLine._start, 
                                _bedLine._end - _bedLine._start, 
                                _step);
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
          _phyloP.processSequence(refSequence, 
                                  _bedLine._blocks[i]._start,
                                  _bedLine._blocks[i]._length,
                                  _step);
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
