/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBEDLINE_H
#define _HALBEDLINE_H

#include <vector>
#include <string>
#include <ostream>
#include "hal.h"

namespace hal {

/** store a block of a BED file (entry from columns 11/12) */
struct BedBlock
{
   hal_index_t _start;
   hal_index_t _length;
   bool operator<(const BedBlock& other) const;
};

/** store a line of a BED file 
 * https://genome.ucsc.edu/FAQ/FAQformat.html#format1
 */
struct BedLine
{
   BedLine();
   virtual ~BedLine();
   std::istream& read(std::istream& is, int version, std::string& lineBuffer);
   std::ostream& write(std::ostream& os, int version=-1);

   std::string _chrName;
   hal_index_t _start;
   hal_index_t _end;
   std::string _name;
   hal_index_t _score;
   char _strand;
   hal_index_t _thickStart;
   hal_index_t _thickEnd;
   hal_index_t _itemR;
   hal_index_t _itemG;
   hal_index_t _itemB;
   std::vector<BedBlock> _blocks;
   std::vector<std::string> _extra;
   int _version;

   // not part of output, but needed to preserve ordering
   hal_index_t _srcStart;
   char _srcStrand;
};

struct BedLineLess
{
  bool operator()(const BedLine& b1, const BedLine& b2) const;
};

struct BedLinePLess
{
  bool operator()(const BedLine* b1, const BedLine* b2) const;
};

struct BedLineSrcLess
{
  bool operator()(const BedLine& b1, const BedLine& b2) const;
};

inline bool BedLinePLess::operator()(const BedLine* b1, const BedLine* b2) const
{
  assert(b1 != NULL && b2 != NULL);
  return BedLineLess().operator()(*b1, *b2);
}

}
#endif
