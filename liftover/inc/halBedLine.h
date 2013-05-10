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
};

/** store a line of a BED file 
 * https://genome.ucsc.edu/FAQ/FAQformat.html#format1
 */
struct BedLine
{
   BedLine();
   virtual ~BedLine();
   void read(std::istream& is, int version, std::string& lineBuffer);
   void write(std::ostream& os, int version);

   std::string _chrName;
   hal_index_t _start;
   hal_index_t _end;
   std::string _name;
   std::string _score;
   char _strand;
   std::string _thickStart;
   std::string _thickEnd;
   std::string _itemRGB;
   std::vector<BedBlock> _blocks;
   std::vector<std::string> _extra;
   int _version;
};

}

#endif
