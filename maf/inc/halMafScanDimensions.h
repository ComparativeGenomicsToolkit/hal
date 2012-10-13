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

   typedef std::pair<hal_size_t, hal_size_t> Interval;
   typedef std::vector<Interval> IntervalList;
   struct Record 
   {
      hal_size_t _length;
      hal_size_t _segments;
      IntervalList _intervals;
   };
   typedef std::map<std::string, Record*> DimMap;

public:

   MafScanDimensions();
   ~MafScanDimensions();
   void scan(const std::string& mafPath);
   const DimMap& getDimensions() const;
   
protected:
   void aLine();
   void sLine();
   void end();
   void updateDimensionsFromBlock();
   void updateDimensionsGlobal();

protected:
      
   DimMap _dimMap;
};

}

#endif
