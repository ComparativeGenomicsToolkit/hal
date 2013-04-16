/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <map>
#include <sstream>
#include "hal.h"
#include "halChain.h"
#include "halBlockViz.h"
#include "halBlockInterpolate.h"

using namespace std;
using namespace hal;


void hal::blockInterpolateSequence(const Sequence* sequence, 
                                   const Genome* target,
                                   hal_size_t start, hal_size_t length,
                                   hal_size_t step)
{
  if (length == 0)
  {
    length = sequence->getSequenceLength() - start;
  }

  set<const Genome*> tgtSet;
  tgtSet.insert(target);
  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(&tgtSet, 0, start,
                                                             start + length-1,
                                                             true);

  hal_index_t absStart = start + sequence->getStartPosition();
  hal_index_t absLast = start + length + sequence->getStartPosition() - 1;
  
  hal_size_t count = 0;
  for (hal_size_t i = 0; i < length; i+= step)
  {
    colIt->toSite(absStart + i, absLast); 
    /* colIt->toRight();
    colIt->toRight();
    colIt->toRight();
    colIt->toRight();
    colIt->toRight();*/
    ++count;
  }
  cout << sequence->getName() << ": " << count << "\n";
}

void hal::blockInterpolateGenome(const Genome* genome, const Sequence* sequence,
                                 const Genome* target,
                                 hal_size_t start, hal_size_t length,
                                 hal_size_t step)
{
  if (sequence != NULL)
  {
    blockInterpolateSequence(sequence, target, start, length, step);
  }
  else
  {
    if (start + length > genome->getSequenceLength())
    {
      stringstream ss;
      ss << "Specified range [" << start << "," << length << "] is"
         << "out of range for genome " << genome->getName() 
         << ", which has length " << genome->getSequenceLength();
      throw (hal_exception(ss.str()));
    }
    if (length == 0)
    {
      length = genome->getSequenceLength() - start;
    }

    SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
    SequenceIteratorConstPtr seqEndIt = genome->getSequenceEndIterator();
    hal_size_t runningLength = 0;
    for (; seqIt != seqEndIt; seqIt->toNext())
    {
      const Sequence* sequence = seqIt->getSequence();
      hal_size_t seqLen = sequence->getSequenceLength();
      hal_size_t seqStart = (hal_size_t)sequence->getStartPosition();

      if (start + length >= seqStart && 
          start < seqStart + seqLen &&
          runningLength < length)
      {
        hal_size_t readStart = seqStart >= start ? 0 : seqStart - start;
        hal_size_t readLen = std::min(seqLen - start, length - runningLength);

        blockInterpolateSequence(sequence, target, readStart, readLen, step);
        runningLength += readLen;
      }
    }
  }
}

