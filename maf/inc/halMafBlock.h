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
#include <map>
#include <string>
#include "hal.h"

namespace hal {

struct MafBlockEntry
{
   std::string _name;
   hal_index_t _start;
   hal_index_t _length;
   char _strand;
   hal_index_t _srcLength;
   std::string _sequence;
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
                  DNAIteratorConstPtr dna);
   void updateEntry(MafBlockEntry* entry, const Sequence* sequence,
                    DNAIteratorConstPtr dna);

   typedef std::multimap<const Sequence*, MafBlockEntry*> Entries;
   Entries _entries;
   Entries::const_iterator _reference;

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

}



#endif
