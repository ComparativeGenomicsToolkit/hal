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
  _topSegment = TopSegmentIteratorPtr();
  _paraTop = TopSegmentIteratorPtr();
  _bottomSegment = BottomSegmentIteratorPtr();
  _refBottom = BottomSegmentIteratorPtr();
  _childIdxMap.clear();
  createGenomes();  
  MafScanner::scan(mafPath, targets);
  initEmptySegments();
  updateRefParseInfo();
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
  _refGenome = _alignment->openGenome(_refName);
  bool newRef = _refGenome == NULL;
  if (newRef)
  {
    if (_alignment->getNumGenomes() != 0)
    {
      throw hal_exception("Cannot add new reference to non-empty alignment");
    }    
    _refGenome = _alignment->addRootGenome(_refName);
  }
  assert(_refGenome != NULL);
  
  // genomes comprise of ranges of sequences in the map.  do a quick
  // scan to add the children of _refGenome.
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
        Sequence::UpdateInfo(sequenceName(i->first), i->second->_numSegments));
    }
    else
    {
      genomeDimensions.push_back(
        Sequence::Info(sequenceName(i->first), i->second->_length, 0, 
                       i->second->_numSegments));
    }    
  }
  if (newRef == false)
  {
    _refGenome->updateBottomDimensions(updateDimensions);
  }
  else
  {
    _refGenome->setDimensions(genomeDimensions);
  }
  if (_refGenome->getNumBottomSegments() > 0)
  {
    _bottomSegment = _refGenome->getBottomSegmentIterator();
    _refBottom = _refGenome->getBottomSegmentIterator();
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
          Sequence::Info(sequenceName(i->first), i->second->_length, 
                         i->second->_numSegments, 0));
      }
      Genome* childGenome = _alignment->openGenome(childName);
      assert(childGenome != NULL);

      childGenome->setDimensions(genomeDimensions);
      if (_topSegment.get() == NULL && childGenome->getNumTopSegments() > 0)
      {
        _topSegment = childGenome->getTopSegmentIterator();
        _paraTop = childGenome->getTopSegmentIterator();
      }
    }
    curRange = getNextSequences(curRange.second);
  }
  
  // update the child idx map
  for (size_t i = 0; i < _refGenome->getNumChildren(); ++i)
  {
    Genome* child = _refGenome->getChild(i);
    _childIdxMap.insert(pair<Genome*, hal_size_t>(child, i));
  }
}

void MafWriteGenomes::convertBlock()
{
  assert(_rows > 0);
  assert(_block[0]._line.length() == _mask.size());
  
  initParaMap();

  for (size_t col = 0; col < _mask.size(); ++col)
  {
    if (_mask[col] == true || col == 0)
    {
      initBlockInfo(col);
      convertSegments(col);
    }
  }
}

void MafWriteGenomes::initBlockInfo(size_t col)
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
      _blockInfo[i]._start = NULL_INDEX;
      _blockInfo[i]._length = 0;
      assert(_dimMap->find(_block[i]._sequenceName) != _dimMap->end());
      _blockInfo[i]._record = _dimMap->find(_block[i]._sequenceName)->second;
      _blockInfo[i]._genome = 
         _alignment->openGenome(genomeName(_block[i]._sequenceName));
      assert(_blockInfo[i]._genome != NULL);
      _blockInfo[i]._skip = false;
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

  _refRow = NULL_INDEX;

  for (size_t i = 0; i < _rows; ++i)
  {
    std::string& line = _block[i]._line;
    Row& row = _block[i];
    RowInfo& rowInfo = _blockInfo[i];
    if (line[col] == '-')
    {
      rowInfo._gaps += last - col;
      _blockInfo[i]._start = NULL_INDEX;
      _blockInfo[i]._length = 0;
    }
    else
    {
      rowInfo._start = row._startPosition + col - rowInfo._gaps;
      rowInfo._length = last - col;
      const StartMap& startMap = rowInfo._record->_startMap;
      const PosSet& posSet = rowInfo._record->_badPosSet;
      if (rowInfo._arrayIndex == NULL_INDEX && rowInfo._skip == false)
      {
        assert(rowInfo._start == row._startPosition);
        assert(rowInfo._gaps <= col);
        StartMap::const_iterator mapIt = startMap.find(rowInfo._start);

        if (mapIt != startMap.end() &&
            mapIt->second._written == 0 &&
            mapIt->second._empty == 0 && 
            posSet.find(FilePosition(_mafFile.tellg(), i)) == posSet.end())
        {
          rowInfo._arrayIndex = mapIt->second._index;
          mapIt->second._written = 1;
          assert(rowInfo._arrayIndex >= 0);
        }
        else
        {
          rowInfo._skip = true;
          if (mapIt->second._empty == 0)
          {
            if (i == (hal_size_t)_refRow)
            {
              stringstream ss;
              ss << "duplication detected in reference at row with start "
                 "position (in + coordinates) " << row._startPosition;
              throw hal_exception(ss.str());
            }
          }
        }
      }
      else if (rowInfo._skip == false)
      {
        ++rowInfo._arrayIndex;
      }
      if (!rowInfo._skip && rowInfo._genome == _refGenome && 
          rowInfo._length > 0 && 
          (_refRow == NULL_INDEX || 
           _block[i]._startPosition < _block[_refRow]._startPosition))
      {
        _refRow = i;
      }
    }
  }
}

void MafWriteGenomes::initParaMap()
{
  for (ParaMap::iterator pIt = _paraMap.begin(); pIt != _paraMap.end(); ++pIt)
  {
    pIt->second.clear();
  }

  for (size_t i = 0; i < _rows; ++i)
  {
    Row& row = _block[i];
    Genome* genome = _alignment->openGenome(genomeName(row._sequenceName));
    assert(genome != NULL);
    Sequence* sequence = genome->getSequence(sequenceName(row._sequenceName));
    assert(sequence != NULL);
    Paralogy para = {sequence->getStartPosition() + row._startPosition, i};
    pair<ParaMap::iterator, bool> res = _paraMap.insert(
      pair<Genome*, ParaSet>(genome, ParaSet()));
    res.first->second.insert(para);    
  }
}

void MafWriteGenomes::convertSegments(size_t col)
{
  // do the reference first
  Sequence* seq;
  size_t numChildren = _refGenome->getNumChildren();
  if (_refRow != NULL_INDEX)
  {
    RowInfo& rowInfo = _blockInfo[_refRow];
    Row& row = _block[_refRow];
    _refBottom->setArrayIndex(_refGenome, rowInfo._arrayIndex);
    seq =  _refGenome->getSequence(sequenceName(row._sequenceName));
    assert(seq != NULL);
    _refBottom->setCoordinates(seq->getStartPosition() + rowInfo._start, 
                               rowInfo._length);
    _refBottom->setTopParseIndex(NULL_INDEX);
    for (size_t i = 0; i < numChildren; ++i)
    {
      _refBottom->setChildIndex(i, NULL_INDEX);
      _refBottom->setChildReversed(i, false);
    }    
    seq->setSubString(
      row._line.substr(col, rowInfo._length), rowInfo._start, rowInfo._length);
  }

  hal_size_t childIndex;
  for (size_t i = 0; i < _rows; ++i)
  {
    if ((hal_index_t)i != _refRow && _blockInfo[i]._length > 0 &&
      _blockInfo[i]._skip == false)
    {
      RowInfo& rowInfo = _blockInfo[i];
      Row& row = _block[i];
      Genome* genome = rowInfo._genome;
      seq = genome->getSequence(sequenceName(row._sequenceName));
      assert(seq != NULL);
      if (genome == _refGenome)
      {
        _bottomSegment->setArrayIndex(rowInfo._genome, rowInfo._arrayIndex);
        _bottomSegment->setCoordinates(seq->getStartPosition() + rowInfo._start,
                                       rowInfo._length);
        _bottomSegment->setTopParseIndex(NULL_INDEX);
        for (size_t i = 0; i < numChildren; ++i)
        {
          _bottomSegment->setChildIndex(i, NULL_INDEX);
          _bottomSegment->setChildReversed(i, false);
        }
      }
      else
      {
        _topSegment->setArrayIndex(rowInfo._genome, rowInfo._arrayIndex);
        _topSegment->setCoordinates(seq->getStartPosition() + rowInfo._start,
                                    rowInfo._length);   

        _topSegment->setNextParalogyIndex(NULL_INDEX);
        _topSegment->setBottomParseIndex(NULL_INDEX);        

        if (_refRow != NULL_INDEX)
        {
          assert(_topSegment->getLength() == _refBottom->getLength());
          childIndex = _childIdxMap.find(rowInfo._genome)->second;
          bool reversed = row._strand != _block[_refRow]._strand;
          _refBottom->setChildIndex(childIndex, rowInfo._arrayIndex);
          _refBottom->setChildReversed(childIndex, reversed);
          _topSegment->setParentIndex(_refBottom->getArrayIndex());
          _topSegment->setParentReversed(reversed);
          updateParalogy(i);
        }
        else
        {
          _topSegment->setParentIndex(NULL_INDEX);
          _topSegment->setParentReversed(false);
        }
      }
      seq->setSubString(
        row._line.substr(col, rowInfo._length), rowInfo._start, 
        rowInfo._length);
    }
  }
}
 
// _top segment needs to be coherent -- we expect its alreay set. 
// only update paralogy if it's already been created (ie previous row)
void MafWriteGenomes::updateParalogy(size_t i)
{
  RowInfo& rowInfo = _blockInfo[i];
  Row& row = _block[i];

  ParaMap::iterator pIt = _paraMap.find(rowInfo._genome);
  assert(pIt != _paraMap.end());
  ParaSet& paraSet = pIt->second;
  if (paraSet.size() > 1)
  {
    Sequence* sequence = rowInfo._genome->getSequence(
      sequenceName(row._sequenceName));
    assert(sequence != NULL);
    Paralogy query = {sequence->getStartPosition() + row._startPosition, 0};
    ParaSet::iterator sIt = paraSet.find(query);
    assert(sIt != paraSet.end());

    ParaSet::iterator next = circularNext(i, paraSet, sIt);
    if (next->_row < i)
    {
      RowInfo& nextInfo = _blockInfo[next->_row];
      assert(nextInfo._length > 0);
      _topSegment->setNextParalogyIndex(nextInfo._arrayIndex);
    }

    ParaSet::iterator prev = circularPrev(i, paraSet, sIt);
    if (prev->_row < i)
    {
      RowInfo& prevInfo = _blockInfo[prev->_row];
      assert(prevInfo._length > 0);
      _paraTop->setArrayIndex(prevInfo._genome, prevInfo._arrayIndex);
      _paraTop->setNextParalogyIndex(rowInfo._arrayIndex);
    }
  }
}

void MafWriteGenomes::initEmptySegments()
{
  for (DimMap::const_iterator dmIt = _dimMap->begin(); dmIt != _dimMap->end();
       ++dmIt)
  {
    Genome* genome = _alignment->openGenome(genomeName(dmIt->first));
    size_t numChildren = genome->getNumChildren();
    Sequence* sequence = genome->getSequence(sequenceName(dmIt->first));
    assert(genome != NULL && sequence != NULL);    
    const StartMap& startMap = dmIt->second->_startMap;

    for (StartMap::const_iterator smIt = startMap.begin();
         smIt != startMap.end(); ++smIt)
    {
      if (smIt->second._written == 0 && smIt->second._count > 0)
      {
        assert(smIt->second._empty == 1);
        const MafScanDimensions::ArrayInfo& arrayInfo = smIt->second;
        hal_size_t startPosition = smIt->first;
        hal_size_t length = 0;
        if ((genome == _refGenome && 
            arrayInfo._index == sequence->getBottomSegmentArrayIndex() +
             sequence->getNumBottomSegments() - 1) ||
            (genome != _refGenome && 
             arrayInfo._index == sequence->getTopSegmentArrayIndex() +
            sequence->getNumTopSegments() - 1))

        {
          length = sequence->getSequenceLength() - startPosition;
        }
        else
        {
          StartMap::const_iterator nextIt = smIt;
          do {
            ++nextIt;
          } while (nextIt != startMap.end() && nextIt->second._count == 0);
          assert(nextIt != startMap.end());
          length = nextIt->first - startPosition;
        }
        if (genome == _refGenome)
        {
          _bottomSegment->setArrayIndex(genome, arrayInfo._index);
          _bottomSegment->setCoordinates(sequence->getStartPosition() + 
                                         startPosition, length);
          _bottomSegment->setTopParseIndex(NULL_INDEX);
          for (size_t i = 0; i < numChildren; ++i)
          {
            _bottomSegment->setChildIndex(i, NULL_INDEX);
            _bottomSegment->setChildReversed(i, false);
          }
        }
        else
        {
          _topSegment->setArrayIndex(genome, arrayInfo._index);
          _topSegment->setCoordinates(sequence->getStartPosition() + 
                                      startPosition, length);
          _topSegment->setBottomParseIndex(NULL_INDEX);
          _topSegment->setParentIndex(NULL_INDEX);
          _topSegment->setParentReversed(false);
          _topSegment->setNextParalogyIndex(NULL_INDEX);
        }
        DNAIteratorPtr dna = sequence->getDNAIterator(startPosition);
        for (hal_size_t i = 0; i < length; ++i)
        {
          dna->setChar('N');
          dna->toRight();
        }
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
    if (_numBlocks % 1000000 == 0)
    {
      cout << "block " << _numBlocks << endl;
      system("date");
    }
  }
}

void MafWriteGenomes::sLine()
{
  MafScanner::Row& row = _block[_rows-1];
  
  // CONVERT TO FORWARD COORDINATES
  if (row._strand == '-')
  {
    row._startPosition =
       row._srcLength - 1 - (row._startPosition + row._length - 1);
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

MafWriteGenomes::ParaSet::iterator 
MafWriteGenomes::circularNext(size_t row, ParaSet& paraSet, ParaSet::iterator i)
{
  ParaSet::iterator j = i;
  do
  {
    ++j;
    if (j == paraSet.end())
    {
      j = paraSet.begin();
    }
    if (_blockInfo[j->_row]._length > 0)
    {
      break;
    }
  } while (i != j);
  return j;
}

MafWriteGenomes::ParaSet::iterator
MafWriteGenomes::circularPrev(size_t row, ParaSet& paraSet, ParaSet::iterator i)
{
  ParaSet::iterator j = i;
  do
  {
    if (j == paraSet.begin())
    {
      j = paraSet.end();
    }
    --j;
    if (_blockInfo[j->_row]._length > 0)
    {
      break;
    }
  } while (i != j);
  return j;
}

void MafWriteGenomes::updateRefParseInfo()
{
  if (_refGenome->getParent() == NULL || _refGenome->getNumChildren() == 0)
  {
    return;
  }
  
  // copied from CactusHalConverter::updateRootParseInfo() in
  // cactus2hal/src/cactusHalConverter.cpp 
  BottomSegmentIteratorPtr bottomIterator = 
     _refGenome->getBottomSegmentIterator();
  TopSegmentIteratorPtr topIterator = _refGenome->getTopSegmentIterator();
  BottomSegmentIteratorConstPtr bend = _refGenome->getBottomSegmentEndIterator();
  TopSegmentIteratorConstPtr tend = _refGenome->getTopSegmentEndIterator();

  while (bottomIterator != bend && topIterator != tend)
  {
    bool bright = false;
    bool tright = false;
    BottomSegment* bseg = bottomIterator->getBottomSegment();
    TopSegment* tseg = topIterator->getTopSegment();
    hal_index_t bstart = bseg->getStartPosition();
    hal_index_t bend = bstart + (hal_index_t)bseg->getLength();
    hal_index_t tstart = tseg->getStartPosition();
    hal_index_t tend = tstart + (hal_index_t)tseg->getLength();
    
    if (bstart >= tstart && bstart < tend)
    {
      bseg->setTopParseIndex(tseg->getArrayIndex());
    }
    if (bend <= tend || bstart == bend)
    {
      bright = true;
    }
        
    if (tstart >= bstart && tstart < bend)
    {
      tseg->setBottomParseIndex(bseg->getArrayIndex());
    }
    if (tend <= bend || tstart == tend)
    {
      tright = true;
    }

    assert(bright || tright);
    if (bright == true)
    {
      bottomIterator->toRight();
    }
    if (tright == true)
    {
      topIterator->toRight();
    }
  }
}
