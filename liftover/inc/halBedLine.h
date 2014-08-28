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

/** hack to store extra fields required to write PSL.  option was never
 * considered in original design and dont have time to refactor properly
 */
struct PSLInfo 
{
   hal_size_t _matches;
   hal_size_t _misMatches;
   hal_size_t _repMatches;
   hal_size_t _nCount;
   hal_size_t _qNumInsert;
   hal_size_t _qBaseInsert;
   hal_size_t _tNumInsert;
   hal_size_t _tBaseInsert;
   std::string _qSeqName;
   hal_size_t _qSeqSize;
   char _qStrand;
   hal_size_t _qEnd;
   hal_size_t _qChromOffset;
   hal_size_t _tSeqSize;
   // absolute sequence coordinates
   std::vector<hal_index_t> _qBlockStarts;
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
   std::ostream& writePSL(std::ostream& os, bool prefixWithName=false);
   bool validatePSL() const;

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

   // hack for optional records that are (in addition to above) required
   // to write psl output.  put in vector so they dont get stored or 
   // copied around if not in use.  vector should never have length  > 1
   std::vector<PSLInfo> _psl;
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
