/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBEDSCANNER_H
#define _HALBEDSCANNER_H

#include <fstream>
#include <string>
#include <deque>
#include <cstdlib>
#include <vector>
#include <string>
#include "hal.h"
#include "halBedLine.h"

namespace hal {

/** Parse a BED file line by line 
 * written independently from the bed export, and it's too much of a 
 * bother to reuse any of that code. */
class BedScanner protected BedLine
{
public:
   BedScanner();
   virtual ~BedScanner();
   virtual void scan(const std::string& bedPath);
   halSize_t getBedVersion() const;
   halSize_t getNumExtraColumns() const;

protected:
   virtual void bedLine() = 0;
   
   std::ifstream _bedFile;

   std::string _chrName;
   hal_index_t _start;
   hal_index_t _end;
   
   
   Block _block;
   size_t _rows;
   Mask _mask;
   hal_size_t _numBlocks;
};

}

#endif
