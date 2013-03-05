/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFSCANNER_H
#define _HALMAFSCANNER_H

#include <fstream>
#include <string>
#include <deque>
#include <cstdlib>
#include <vector>
#include <string>
#include "hal.h"

namespace hal {

/** Parse a MAF file line by line 
 * written independently from the maf export, and it's too much of a 
 * bother to reuse any of that code. */
class MafScanner
{
public:
   MafScanner();
   virtual ~MafScanner();
   virtual void scan(const std::string& mafPath, 
                     const std::set<std::string>& targetSet);
   hal_size_t getNumBlocks() const { return _numBlocks; }
   static std::string genomeName(const std::string& fullName);
   static std::string sequenceName(const std::string& fullName);

   struct Row {
      std::string _sequenceName;
      hal_size_t _startPosition;
      hal_size_t _length;
      char _strand;
      hal_size_t _srcLength;
      std::string _line;
   };
   typedef std::vector<Row> Block;
   typedef std::vector<bool> Mask;

protected:
   virtual void aLine() = 0;
   virtual void sLine() = 0;
   virtual void end() = 0;
   void nextLine();
   void updateMask();


   std::ifstream _mafFile;
   std::set<std::string> _targets;
   
   Block _block;
   size_t _rows;
   Mask _mask;
   hal_size_t _numBlocks;
};

}

#endif
