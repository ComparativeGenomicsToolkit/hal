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
#include "sonLib.h"

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
   // add this because _sequence is no longer assumed to 
   // be unique
   const Genome* _genome;
   // The node corresponding to this entry (if we are printing trees)
   stTree *_tree;
};

class MafBlock
{
public:
   static const hal_index_t defaultMaxLength;

   MafBlock(hal_index_t maxLength = defaultMaxLength);
   ~MafBlock();

   void initBlock(ColumnIteratorConstPtr col, bool fullNames, bool printTree);
   void appendColumn(ColumnIteratorConstPtr col);
   bool canAppendColumn(hal::ColumnIteratorConstPtr col);
   void setMaxLength(hal_index_t maxLen);
   
protected:
   
   void resetEntries();
   void initEntry(MafBlockEntry* entry, const Sequence* sequence,
                  DNAIteratorConstPtr dna, bool clearSequence = true);
   void updateEntry(MafBlockEntry* entry, const Sequence* sequence,
                    DNAIteratorConstPtr dna);
   std::string getName(const Sequence* sequence) const;
   stTree *buildTree(ColumnIteratorConstPtr colIt, bool modifyEntries);
   void buildTreeR(BottomSegmentIteratorConstPtr botIt, stTree *tree, bool modifyEntries);
   stTree *getTreeNode(SegmentIteratorConstPtr segIt, bool modifyEntries);

   std::ostream& printBlock(std::ostream& os) const;
   std::ostream& printBlockWithTree(std::ostream& os) const;

   typedef std::multimap<const Sequence*, MafBlockEntry*, 
                         ColumnIterator::SequenceLess> Entries;
   Entries _entries;
   Entries::const_iterator _reference;
   std::vector<MafBlockString*> _stringBuffers;
   hal_index_t _maxLength;
   hal_index_t _refIndex;
   bool _fullNames;
   bool _printTree;
   stTree *_tree;

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
  _buffers(buffers), _lastUsed(0), _genome(NULL)
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

inline std::string MafBlock::getName(const Sequence* sequence) const
{
  return _fullNames ? sequence->getFullName() : sequence->getName();
}

inline void MafBlock::setMaxLength(hal_index_t maxLen) 
{
  _maxLength = maxLen;
}

}



#endif
