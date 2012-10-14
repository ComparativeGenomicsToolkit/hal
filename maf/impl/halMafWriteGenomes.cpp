/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "hal.h"
#include "halMafWriteGenomes.h"

using namespace std;
using namespace hal;


MafWriteGenomes::MafWriteGenomes() : MafScanner()
{

}

MafWriteGenomes::~MafWriteGenomes()
{

}

void MafWriteGenomes::convert(const string& mafPath,
                              const string& refGenomeName,
                              const set<string>& targets,
                              const DimMap& dimMap,
                              AlignmentPtr alignment)
{
  _refName = refGenomeName;
  _dimMap = &dimMap;
  _alignment = alignment;

  createGenomes();  
  MafScanner::scan(mafPath, targets);
}

MafWriteGenomes::MapRange MafWriteGenomes::getRefSequences() const
{
  DimMap::const_iterator i = _dimMap->lower_bound(_refName);
  for (; i != _dimMap->end(); ++i)
  {
    if (genomeName(i->first) == _refName)
    {
      break;
    }
  }
  
  if (i == _dimMap->end())
  {
    stringstream ss;
    ss << "Reference genome " << _refName << " was not found in the MAF file.  "
       "Ie, no sequence in name in the form " << _refName << ".something was "
       "found.";
    throw hal_exception(ss.str());
  }

  DimMap::const_iterator j = i;
  for (; j != _dimMap->end(); ++j)
  {
    if (genomeName(j->first) != _refName)
    {
      break;
    }
  }
  
  return MapRange(i, j);
}

MafWriteGenomes::MapRange MafWriteGenomes::getNextSequences(
  DimMap::const_iterator jprev) const
{
  DimMap::const_iterator j = jprev;
  if (j == _dimMap->end())
  {
    return MapRange(j, j);
  }
  string name = genomeName(j->first);
  for (; j != _dimMap->end(); ++j)
  {
    if (genomeName(j->first) != name)
    {
      break;
    }
  }
  return MapRange(jprev, j);
}

void MafWriteGenomes::createGenomes()
{
  // need to create the tree before we do anything
  Genome* refGenome = _alignment->openGenome(_refName);
  bool newRef = refGenome == NULL;
  if (newRef)
  {
    if (_alignment->getNumGenomes() != 0)
    {
      throw hal_exception("Cannot add new reference to non-empty alignment");
    }    
    refGenome = _alignment->addRootGenome(_refName);
  }
  assert(refGenome != NULL);
  
  // genomes comprise of ranges of sequences in the map.  do a quick
  // scan to add the children of refGenome.
  MapRange refRange = getRefSequences();
  MapRange curRange = getNextSequences(_dimMap->begin());
  while (curRange.first != _dimMap->end())
  {
    if (curRange != refRange)
    {
      // note, we do not add a branch length.  will need to add option
      // to maf2hal where a tree can be given as well just for branch lengths
      _alignment->addLeafGenome(genomeName(curRange.first->first), _refName, 1);
    }
    curRange = getNextSequences(curRange.second);
  }

  // do the ref genome dimenions
  vector<Sequence::Info> genomeDimensions;
  vector<Sequence::UpdateInfo> updateDimensions;
  for (DimMap::const_iterator i = refRange.first; i != refRange.second; ++i)
  {
    if (newRef == false)
    {
      updateDimensions.push_back(
        Sequence::UpdateInfo(i->first, i->second->_startMap.size()));
    }
    else
    {
      genomeDimensions.push_back(
        Sequence::Info(i->first, i->second->_length, 0, 
                       i->second->_startMap.size()));
    }    
  }
  if (newRef == false)
  {
    refGenome->updateBottomDimensions(updateDimensions);
  }
  else
  {
    refGenome->setDimensions(genomeDimensions);
  }

  // do the child genome dimensions
  curRange = getNextSequences(_dimMap->begin());
  while (curRange.first != _dimMap->end())
  {
    if (curRange != refRange)
    {
      string childName = genomeName(curRange.first->first);
      genomeDimensions.clear();
      for (DimMap::const_iterator i = curRange.first; i != curRange.second; ++i)
      {
        genomeDimensions.push_back(
          Sequence::Info(i->first, i->second->_length, 
                         i->second->_startMap.size(), 0));
      }
      Genome* childGenome = _alignment->openGenome(childName);
      assert(childGenome != NULL);
      childGenome->setDimensions(genomeDimensions);
    }
    curRange = getNextSequences(curRange.second);
  }
}

void MafWriteGenomes::convertBlock()
{
  assert(_rows > 0);
  assert(_block[0]._line.length() == _mask.size());
  initArrayIndexes(0);

  for (size_t col = 1; col < _mask.size(); ++col)
  {
    if (_mask[col] == true)
    {
      initArrayIndexes(col);
    }
  }
}

void MafWriteGenomes::initArrayIndexes(size_t col)
{
  if (col == 0)
  {
    if (_blockInfo.size() < _rows)
    {
      _blockInfo.resize(_rows);
    }
    for (size_t i = 0; i < _rows; ++i)
    {
      _blockInfo[i]._arrayIndex = NULL_INDEX;
      _blockInfo[i]._gaps = 0;
      assert(_dimMap->find(_block[i]._sequenceName) != _dimMap->end());
      _blockInfo[i]._record = _dimMap->find(_block[i]._sequenceName)->second;
    }
  }
  else
  {
    assert(_blockInfo.size() >= _rows);
  }

  // our range will be [col, last)
  size_t last = col + 1;
  while (last < _mask.size() && _mask[last] == false)
  {
    ++last;
  }

  for (size_t i = 0; i < _rows; ++i)
  {
    std::string& line = _block[i]._line;
    Row& row = _block[i];
    RowInfo& rowInfo = _blockInfo[i];
    if (line[col] == '-')
    {
      rowInfo._gaps += last - col;
    }
    else
    {
      const StartMap& startMap = rowInfo._record->_startMap;
      if (rowInfo._arrayIndex == NULL_INDEX)
      {
        assert(rowInfo._gaps <= col);
        hal_size_t startPos = row._startPosition + col - rowInfo._gaps;
        StartMap::const_iterator mapIt = startMap.find(startPos);
        assert(mapIt != startMap.end());
        rowInfo._arrayIndex = mapIt->second;
      }
      else
      {
        assert(startMap.find(row._startPosition + col - rowInfo._gaps) 
               != startMap.end());
        assert(startMap.find(row._startPosition + col - rowInfo._gaps)->second
               == (hal_size_t)rowInfo._arrayIndex);
      }
    }
  }
}

void MafWriteGenomes::aLine()
{
  assert(_rows <= _block.size());
  if (_rows > 0)
  {
    convertBlock();
  }
}

void MafWriteGenomes::sLine()
{
  MafScanner::Row& row = _block[_rows-1];
  
  // CONVERT TO FORWARD COORDINATES
  if (row._strand == '-')
  {
    row._startPosition = row._srcLength - 1 - row._length;
    reverseComplement(row._line);
  }
}

void MafWriteGenomes::end()
{
  assert(_rows <= _block.size());
  if (_rows > 0)
  {
    convertBlock();
  }
}
