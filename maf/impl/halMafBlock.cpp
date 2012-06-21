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

void MafBlock::setFirstColumn(ColumnIteratorConstPtr col)
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
    }    
  }
  appendColumn(col);
}

void MafBlock::appendColumn(ColumnIteratorConstPtr col)
{
  const ColumnMap* colMap = col->getColumnMap();
  size_t idx = 1;
  bool refFound = false;
  const Sequence* refSequence = col->getReferenceSequence();
  string name;
  MafBlockEntry* entry = NULL;
  // for every sequence
  for (ColumnMap::const_iterator i = colMap->begin(); i != colMap->end(); ++i)
  {
    const Sequence* sequence = i->first;
    name = sequence->getName();
 
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
             (*j)->getArrayIndex() == entry->_length - entry->_start);
      assert((*j)->getReversed() == false || 
             (*j)->getArrayIndex() == entry->_start - entry->_length);
      entry->_sequence.append(1, (*j)->getChar());
    }
    
    // UPDATE GAPPED ENTRIES
    while (_entries[idx]._name == name)
    {
      entry->_sequence.append("-");
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
  string name;

  // for every sequence
  for (ColumnMap::const_iterator i = colMap->begin(); i != colMap->end(); ++i)
  {
    const Sequence* sequence = i->first;
    name = sequence->getName();

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
        entry = &_entries[idx++];
      }
      if (entry->_name != name)
      {
        return false;
      }
      if (entry->_strand == '+')
      {
        if ((*j)->getArrayIndex() - entry->_start != entry->_length)
        {
          return false;
        }
      }
      else
      {
        assert(entry->_strand == '-');
        if (entry->_start - (*j)->getArrayIndex() != entry->_length)
        {
          return false;
        }        
      }
    }

    while (_entries[idx]._name == name)
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
  os << endl;
  return os;
}

istream& hal::operator>>(istream& is, MafBlock& mafBlock)
{
  return is;
}
