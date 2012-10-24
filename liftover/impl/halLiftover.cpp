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
    bool lineRead  = readBedLine();
    if (lineRead == true)
    {
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
}

void Liftover::liftInterval()
{  
  map<const Sequence*, PositionCache*> posCacheMap;
  _colIt = _srcSequence->getColumnIterator(&_tgtSet, 0, _inStart, _inEnd - 1);
  do 
  {
    const ColumnMap* cMap = _colIt->getColumnMap();
    for (ColumnMap::const_iterator i = cMap->begin(); i != cMap->end(); ++i)
    {
      if (i->first->getGenome() == _tgtGenome)
      {
        const DNASet* dSet = i->second;
        const Sequence* seq = i->first;
        for (DNASet::const_iterator j = dSet->begin(); j != dSet->end(); ++j)
        {
          pair<map<const Sequence*, PositionCache*>::iterator, bool> res =
             posCacheMap.insert(pair<const Sequence*, PositionCache*>(
                                  seq, NULL));
          if (res.second == true)
          {
            res.first->second = new PositionCache();
          }
          res.first->second->insert((*j)->getArrayIndex());
        }
      }
    }
    _colIt->toRight();
  } 
  while (_colIt->lastColumn() == false);

  map<const Sequence*, PositionCache*>::iterator pcmIt;
  for (pcmIt = posCacheMap.begin(); pcmIt != posCacheMap.end(); ++pcmIt)
  {
    const Sequence* seq = pcmIt->first;
    hal_size_t seqStart = seq->getStartPosition();
    PositionCache* posCache = pcmIt->second;
    const IntervalSet* iSet = posCache->getIntervalSet();
    _outName = seq->getName();
    for (IntervalSet::const_iterator k = iSet->begin(); k != iSet->end(); ++k)
    {
      _outStart = k->second - seqStart;
      _outEnd = k->first + 1 - seqStart;
      writeBedLine();
    }
    delete posCache;
  }
}

bool Liftover::readBedLine()
{
  _buffer.clear();
  _data.clear();
  getline(*_inputFile, _buffer);
  istringstream ss(_buffer);
  ss >> _inName;
  if (ss.bad() || ss.eof() || _inName[0] == '#')
  {
    return false;
  }
  ss >> _inStart >> _inEnd;
  string token;
  while (!ss.bad() && !ss.eof())
  {
    ss >> token;
    _data.push_back(token);
  }
  return !ss.bad();
}

void Liftover::writeBedLine()
{
  *_outputFile << _outName << '\t' << _outStart << '\t' << _outEnd;
  for (size_t i = 0; i < _data.size(); ++i)
  {
    *_outputFile << '\t' << _data[i];
  }
  *_outputFile << '\n';
}
