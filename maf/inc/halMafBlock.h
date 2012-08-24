/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFBLOCK_H
#define _HALMAFBLOCK_H

#include <iostream>
#include <string>
#include <deque>
#include <cstdlib>
#include <map>
#include <string>
#include "hal.h"

namespace hal {

// used to use std::string but appending was too slow
// though sometimes i think instruments' cpu profiler has been lying
// to me and i did all this for nothing. 
struct MafBlockString
{
   MafBlockString();
   ~MafBlockString();
   void append(char c);
   void clear();
   char* str();
   char* _buf;
   size_t _cap;
   size_t _len;
};

struct MafBlockEntry
{
   // we hack to keep a global buffer list to reduce 
   // allocs and frees as entries get created and destroyed
   MafBlockEntry(std::vector<MafBlockString*>& buffers);
   ~MafBlockEntry();
   
   std::vector<MafBlockString*>& _buffers;
   std::string _name;
   hal_index_t _start;
   hal_index_t _length;
   char _strand;
   short _lastUsed;
   hal_index_t _srcLength;
   MafBlockString* _sequence;
};

class MafBlock
{
public:

   MafBlock();
   ~MafBlock();

   void initBlock(ColumnIteratorConstPtr col);
   void appendColumn(ColumnIteratorConstPtr col);
   bool canAppendColumn(hal::ColumnIteratorConstPtr col);
   
protected:
   
   void resetEntries();
   void initEntry(MafBlockEntry* entry, const Sequence* sequence,
                  DNAIteratorConstPtr dna, bool clearSequence = true);
   void updateEntry(MafBlockEntry* entry, const Sequence* sequence,
                    DNAIteratorConstPtr dna);

   typedef std::multimap<const Sequence*, MafBlockEntry*, 
                         ColumnIterator::SequenceLess> Entries;
   Entries _entries;
   Entries::const_iterator _reference;
   std::vector<MafBlockString*> _stringBuffers;

   typedef hal::ColumnIterator::ColumnMap ColumnMap;
   typedef hal::ColumnIterator::DNASet DNASet;
   friend std::ostream& operator<<(std::ostream& os, 
                                   const hal::MafBlock& mafBlock);
   friend std::istream& operator>>(std::istream& is, hal::MafBlock& mafBlock);
};

std::ostream& operator<<(std::ostream& os, 
                        const hal::MafBlockEntry& mafBlockEntry);
std::istream& operator>>(std::istream& is, hal::MafBlockEntry& mafBlockEntry);
std::ostream& operator<<(std::ostream& os, const hal::MafBlock& mafBlock);
std::istream& operator>>(std::istream& is, hal::MafBlock& mafBlock);

inline MafBlockString::MafBlockString() :
  _buf(NULL),
  _cap(0),
  _len(0)
{
}

inline MafBlockString::~MafBlockString()
{
  free(_buf);
}

inline void MafBlockString::append(char c)
{
  assert(_cap >= _len);
  if (_cap == 0)
  {
    _cap = 127;
    _buf = (char*)malloc((_cap + 1) * sizeof(char));
  }
  else if (_cap == _len)
  {
    _cap = (_cap + 1) * 2 - 1;
    _buf = (char*)realloc(_buf, _cap + 1); 
  }
  _buf[_len++] = c;
}

inline void MafBlockString::clear()
{
  if (_cap > 0)
  {
    _len = 0;
  }
}

inline char* MafBlockString::str()
{
  assert(_cap > 0);
  _buf[_len] = '\0';
  return _buf;
}

inline MafBlockEntry::MafBlockEntry(std::vector<MafBlockString*>& buffers) : 
  _buffers(buffers), _lastUsed(0) 
{
  if (_buffers.empty() == false) 
  {
    _sequence = _buffers.back();
    _buffers.pop_back();
    _sequence->clear();
  }   
  else
  {
    _sequence = new MafBlockString();
  }
}

inline MafBlockEntry::~MafBlockEntry() 
{
  if (_sequence != NULL) 
  {
    // do we need to keep track of a maximum size here?  probably not,
    // since the number of entries is bounded by the number in use
    // at any one time which shouldn't be too outrageous.  still maybe
    // should have a hard limit for sanity. 
    _buffers.push_back(_sequence);
  }
}


}



#endif