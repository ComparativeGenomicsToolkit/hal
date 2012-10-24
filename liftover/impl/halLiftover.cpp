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

Liftover::Liftover() : _inputFile(NULL), _outputFile(NULL),
                       _srcGenome(NULL), _tgtGenome(NULL)
{

}

Liftover::~Liftover()
{

}

void Liftover::convert(AlignmentConstPtr alignment,
                       const Genome* srcGenome,
                       ifstream* inputFile,
                       const Genome* tgtGenome,
                       ofstream* outputFile)
{
  _srcGenome = srcGenome;
  _inputFile = inputFile;
  _tgtGenome = tgtGenome;
  _outputFile = outputFile;
  _missedSet.clear();
  _tgtSet.clear();
  assert(_srcGenome && _inputFile && tgtGenome && outputFile);

  _tgtSet.insert(tgtGenome);

  while (!inputFile->bad() && !inputFile->eof())
  {
    readBedLine();
    _srcSequence = _srcGenome->getSequence(_inName);
    if (_srcSequence == NULL)
    {
      pair<set<string>::iterator, bool> result = _missedSet.insert(_inName);
      if (result.second == true)
      {
        std::cerr << "Unable to find sequence " << _inName << " in genome "
                  << srcGenome->getName() << endl;
      }
    }
    
    else if (_inEnd >= _srcSequence->getSequenceLength())
    {
      std::cerr << "Skipping interval with endpoint " << _inEnd 
                << "because sequence " << _inName << " has length " 
                << _srcSequence->getSequenceLength() << endl;
    }
    else
    {
      liftInterval();
    }
  }
}

void Liftover::liftInterval()
{  
  _colIt = _srcSequence->getColumnIterator(&_tgtSet, 0, _inStart, _inEnd - 1);
  do 
  {
    const ColumnMap* colMap = _colIt->getColumnMap();
    _colIt->toRight();
  } 
  while (_colIt->lastColumn() == false);
  
}

void Liftover::readBedLine()
{

}

void Liftover::writeBedLine()
{

}
