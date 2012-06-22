/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include "halMafBlock.h"

using namespace std;
using namespace hal;

MafBlock::MafBlock()
{

}

MafBlock::~MafBlock()
{

}

void MafBlock::initBlock(ColumnIteratorConstPtr col)
{
  _entries.clear();
  const ColumnMap* colMap = col->getColumnMap();
  _lastSize = colMap->size();
  bool refAdded = false;
  const Sequence* refSequence = col->getReferenceSequence();
  // for every sequence
  for (ColumnMap::const_iterator i = colMap->begin(); i != colMap->end(); ++i)
  {
    const Sequence* sequence = i->first;
    // for every base in each sequence (<=1 unless duplication)
    for (DNASet::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
    {
      MafBlockEntry* entry = NULL;
      if (refAdded == false && sequence == refSequence)
      {
        refAdded = true;
        _entries.push_front(MafBlockEntry());
        entry = &_entries.front();
      }
      else
      {
        _entries.push_back(MafBlockEntry());
        entry = &_entries.back();
      }
      entry->_name = sequence->getName();
      entry->_start = (*j)->getArrayIndex() - sequence->getStartPosition();
      entry->_length = 0;
      entry->_strand = (*j)->getReversed() ? '-' : '+';
      entry->_srcLength = (hal_index_t)sequence->getSequenceLength();
      if ((*j)->getReversed())
      {
        entry->_start = entry->_srcLength - entry->_start;
      }
    }    
  }
}

void MafBlock::appendColumn(ColumnIteratorConstPtr col)
{
  const ColumnMap* colMap = col->getColumnMap();
  size_t idx = 1;
  bool refFound = false;
  const Sequence* refSequence = col->getReferenceSequence();
  MafBlockEntry* entry = NULL;
  // for every sequence
  for (ColumnMap::const_iterator i = colMap->begin(); i != colMap->end(); ++i)
  {
    const Sequence* sequence = i->first;
    hal_index_t sequenceStart = sequence->getStartPosition();
    const string& name = sequence->getName();
 
    // UPDATE NON-GAPPED ENTRIES
    for (DNASet::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
    {
      if (refFound == false && sequence == refSequence)
      {
        refFound = true;
        entry = &_entries.front();
      }
      else
      {
        entry = &_entries[idx++];
      }
      assert(entry->_name == sequence->getName());
      assert(entry->_strand == (*j)->getReversed() ? '-' : '+');
      assert(entry->_srcLength == (hal_index_t)sequence->getSequenceLength());
      ++entry->_length;
      assert((*j)->getReversed() == true || 
             (*j)->getArrayIndex() - sequenceStart == 
             entry->_start + entry->_length - 1);
      assert((*j)->getReversed() == false || 
             entry->_srcLength - (*j)->getArrayIndex() + sequenceStart == 
             entry->_start + entry->_length - 1);
      entry->_sequence.append(1, (*j)->getChar());
    }
    
    // UPDATE GAPPED ENTRIES
    while (idx < _entries.size() && _entries[idx]._name == name)
    {
      _entries[idx]._sequence.append("-");
      ++idx;
    } 
  }
}

// Q: When can we append a column? 
// A: When for every sequence already in the column, the new column
//    has either a gap or a contigugous base.  The new column also has
//    no new sequences.  
bool MafBlock::canAppendColumn(ColumnIteratorConstPtr col)
{
  const ColumnMap* colMap = col->getColumnMap();
  if (colMap->size() != _lastSize)
  {
    return false;
  }
  size_t idx = 1;
  bool refFound = false;
  const Sequence* refSequence = col->getReferenceSequence();
  MafBlockEntry* entry = NULL;

  // for every sequence
  for (ColumnMap::const_iterator i = colMap->begin(); i != colMap->end(); ++i)
  {
    const Sequence* sequence = i->first;
    hal_index_t sequenceStart = sequence->getStartPosition();
    hal_index_t pos = 0;
    const string& name = sequence->getName();

    // for every base in each sequence (<=1 unless duplication)
    for (DNASet::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
    {
      if (refFound == false && sequence == refSequence)
      {
        refFound = true;
        entry = &_entries.front();
      }
      else
      {
        if (idx >= _entries.size())
        {
          return false;
        }
        entry = &_entries[idx++];
      }
      if (entry->_name != name)
      {
        return false;
      }

      // position on forward strand relative to start of sequence
      pos = (*j)->getArrayIndex() - sequenceStart;
      if (entry->_strand == '-')
      {
        assert(entry->_strand == '-');
        // position on reverse strand relative to end of sequence
        pos = entry->_srcLength - pos;
      }

      if (pos - entry->_start != entry->_length)
      {
        return false;
      }
    }

    while (idx < _entries.size() && _entries[idx]._name == name)
    {
      ++idx;
    } 
  } 
  return true;
}

ostream& hal::operator<<(ostream& os, const MafBlockEntry& mafBlockEntry)
{
  os << "s\t" << mafBlockEntry._name << '\t' << mafBlockEntry._start << '\t'
     << mafBlockEntry._length << '\t' << mafBlockEntry._strand << '\t' 
     << mafBlockEntry._srcLength << '\t' << mafBlockEntry._sequence << '\n';
  return os;
}

istream& hal::operator>>(istream& is, MafBlockEntry& mafBlockEntry)
{
  is >> mafBlockEntry._name >> mafBlockEntry._start 
     >> mafBlockEntry._length >> mafBlockEntry._strand 
     >> mafBlockEntry._srcLength >> mafBlockEntry._sequence;
  assert(mafBlockEntry._strand == '+' || mafBlockEntry._strand == '-');
  return is;
}

ostream& hal::operator<<(ostream& os, const MafBlock& mafBlock)
{
  os << "a\n";
  for (size_t i = 0; i < mafBlock._entries.size(); ++i)
  {
    os << mafBlock._entries[i];
  }
  return os;
}

istream& hal::operator>>(istream& is, MafBlock& mafBlock)
{
  return is;
}
