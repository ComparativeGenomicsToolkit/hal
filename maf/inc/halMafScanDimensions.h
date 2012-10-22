/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFSCANDIMENSIONS_H
#define _HALMAFSCANDIMENSIONS_H

#include <iostream>
#include <string>
#include <deque>
#include <cstdlib>
#include <map>
#include <string>
#include "halMafScanner.h"

namespace hal {

/** Parse a MAF file line by line, getting some dimension stats
 * and maybe checking for some errros. */
class MafScanDimensions : public MafScanner
{
public:
   // map start position to array index
   struct ArrayInfo 
   {
      hal_size_t _index : 42;
      hal_size_t _count : 20;
      hal_size_t _empty : 1;
      mutable hal_size_t _written : 1;
   };
   typedef std::map<hal_size_t, ArrayInfo> StartMap;

   typedef std::pair<std::streampos, size_t> FilePosition;
   typedef std::set<FilePosition> PosSet;

   struct Record 
   {
      hal_size_t _length;
      hal_size_t _numSegments;
      StartMap _startMap;
      PosSet _badPosSet;
   };
   typedef std::map<std::string, Record*> DimMap;

public:

   MafScanDimensions();
   ~MafScanDimensions();
   void scan(const std::string& mafPath, 
             const std::set<std::string>& targetSet);
   const DimMap& getDimensions() const;
   
protected:
   void aLine();
   void sLine();
   void end();
   void updateDimensionsFromBlock();
   void updateArrayIndices();

protected:
      
   DimMap _dimMap;
};

}

#endif
