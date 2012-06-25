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
  for (Entries::iterator i = _entries.begin(); i != _entries.end(); ++i)
  {
    delete i->second;
  }
}

void MafBlock::resetEntries()
{
  Entries::iterator i = _entries.begin();
  Entries::iterator next;
  MafBlockEntry* e;
  while (i != _entries.end())
  {
    next = i;
    ++next;
    e = i->second;
    if (e->_start == NULL_INDEX)
    {
      // Never got written to in the last block, so we give up
      // and remove it.
      delete i->second;
      _entries.erase(i);
    }
    else
    {
      assert (e->_length > 0);
      assert (e->_name == i->first->getName());
      // Rest block information but leave sequence information so we
      // can reuse it. 
      e->_start = NULL_INDEX;
      e->_strand = '+';
      e->_length = 0;
      e->_sequence.clear();
    }
    i = next;
  }  
}

void MafBlock::initEntry(MafBlockEntry* entry, const Sequence* sequence, 
                         DNAIteratorConstPtr dna)
{
  if (entry->_name != sequence->getName())
  {
    // replace genearl sequence information
    entry->_name = sequence->getName();
    entry->_srcLength = (hal_index_t)sequence->getSequenceLength();
  }
  if (dna.get())
  {
    // update start position from the iterator
    entry->_start = dna->getArrayIndex() - sequence->getStartPosition();
    entry->_length = 0;
    entry->_strand = dna->getReversed() ? '-' : '+';

    if (dna->getReversed())
    {
      entry->_start = entry->_srcLength - entry->_start;
    }
  }
  else
  {
    // no start position, so we wait till next time
    entry->_start = NULL_INDEX;
    entry->_length = 0;
    entry->_strand = '+';
  }
  entry->_sequence.clear();
}

void MafBlock::updateEntry(MafBlockEntry* entry, const Sequence* sequence,
                           DNAIteratorConstPtr dna)
{
  if (dna.get() != NULL)
  {
    if (entry->_start == NULL_INDEX)
    {
      string prep = entry->_sequence;
      initEntry(entry, sequence, dna);
      entry->_sequence.swap(prep);
      entry->_sequence.append(prep);
    }
    assert(entry->_name == sequence->getName());
    assert(entry->_strand == dna->getReversed() ? '-' : '+');
    assert(entry->_srcLength == (hal_index_t)sequence->getSequenceLength());
    
    ++entry->_length;
    
    assert(dna->getReversed() == true || 
           dna->getArrayIndex() - sequence->getStartPosition() == 
           entry->_start + entry->_length - 1);
    assert(dna->getReversed() == false || 
           entry->_srcLength - dna->getArrayIndex() + 
           sequence->getStartPosition()  == 
           entry->_start + entry->_length - 1);

    entry->_sequence.push_back(dna->getChar());
  }
  else
  {
    entry->_sequence.append("-");
  }
}

void MafBlock::initBlock(ColumnIteratorConstPtr col)
{
  resetEntries();
  const ColumnMap* colMap = col->getColumnMap();
  Entries::iterator e = _entries.begin();
  ColumnMap::const_iterator c = colMap->begin();
  DNASet::const_iterator d;
  const Sequence* sequence;
  
  for (; c != colMap->end(); ++c)
  {
    sequence = c->first;

    // No DNA Iterator for this sequence.  We just give it an empty
    // entry
    if (c->second->empty())
    {
      e = _entries.lower_bound(sequence);
      if (e == _entries.end() || e->first != sequence)
      {
        MafBlockEntry* entry = new MafBlockEntry();
        initEntry(entry, sequence, DNAIteratorConstPtr());
        e = _entries.insert(Entries::value_type(sequence, entry));
      }
      else
      {
        assert (e->first == sequence);
        assert (e->second->_name == sequence->getName());
        initEntry(e->second, sequence, DNAIteratorConstPtr());
      }
    }

    else
    {
      for (d = c->second->begin(); d != c->second->end(); ++d)
      {
        // search for c's sequence in _entries.  
        // we conly call find() once.  afterwards we just move forward
        // in the map since they are both sorted by the same key. 
        if (e == _entries.begin())
        {
          e = _entries.lower_bound(sequence);
          if (e == _entries.end() || e->first != sequence)
          {
            e = _entries.end();
          }
        }
        else
        {
          ++e;
          while (e->first != c->first && e != _entries.end())
          {
            ++e;
          }
        }
        if (e == _entries.end())
        {
          MafBlockEntry* entry = new MafBlockEntry();
          initEntry(entry, sequence, *d);
          assert(entry->_name == sequence->getName());
          e = _entries.insert(Entries::value_type(sequence, entry));
        }
        else
        {
          initEntry(e->second, sequence, *d);
        }
      }
    }
  }
}    

void MafBlock::appendColumn(ColumnIteratorConstPtr col)
{
  const ColumnMap* colMap = col->getColumnMap();
  Entries::iterator e = _entries.begin();
  ColumnMap::const_iterator c = colMap->begin();
  DNASet::const_iterator d;
  const Sequence* sequence;

  for (; c != colMap->end(); ++c)
  {
    sequence = c->first;
    for (d = c->second->begin(); d != c->second->end(); ++d)
    {
      while (e->first != sequence && e != _entries.end())
      {
        updateEntry(e->second, NULL, DNAIteratorConstPtr());
        ++e;
      }
      assert(e != _entries.end());
      assert(e->first == sequence);
      assert(e->second->_name == sequence->getName());
      updateEntry(e->second, sequence, *d);
      ++e;
    }
  }
  
  for (; e != _entries.end(); ++e)
  {
    updateEntry(e->second, NULL, DNAIteratorConstPtr());
  }
}


// Q: When can we append a column? 
// A: When for every sequence already in the column, the new column
//    has either a gap or a contigugous base.  The new column also has
//    no new sequences.  
bool MafBlock::canAppendColumn(ColumnIteratorConstPtr col)
{
  const ColumnMap* colMap = col->getColumnMap();
  Entries::iterator e = _entries.begin();
  ColumnMap::const_iterator c;
  DNASet::const_iterator d;
  const Sequence* sequence;
  MafBlockEntry* entry;
  hal_index_t sequenceStart;
  hal_index_t pos;
 
  for (c = colMap->begin(); c != colMap->end(); ++c)
  {
    sequence = c->first;
    sequenceStart = sequence->getStartPosition();

    for (d = c->second->begin(); d != c->second->end(); ++d)
    {
      while (e->first != sequence && e != _entries.end())
      {
        ++e;
      }
      if (e == _entries.end())
      {
        return false;
      }
      else
      {
        entry = e->second;
        assert(e->first == sequence);
        assert(entry->_name == sequence->getName());
        if (entry->_start != NULL_INDEX)
        {
          pos = (*d)->getArrayIndex() - sequenceStart;
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
        ++e;
      }
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

// todo: fast way of reference first. 
ostream& hal::operator<<(ostream& os, const MafBlock& mafBlock)
{
  os << "a\n";
  for (MafBlock::Entries::const_iterator e = mafBlock._entries.begin();
       e != mafBlock._entries.end(); ++e)
  {
    if (e->second->_start != NULL_INDEX)
    {
      os << *e->second;
    }
  }
  return os;
}

istream& hal::operator>>(istream& is, MafBlock& mafBlock)
{
  return is;
}
