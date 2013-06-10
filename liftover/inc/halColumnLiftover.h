/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNLIFTOVER_H
#define _HALCOLUMNLIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "halLiftover.h"

namespace hal {

class ColumnLiftover : public Liftover
{
public:
   
   ColumnLiftover();
   virtual ~ColumnLiftover();
                   
protected:

   void liftInterval(BedList& mappedBedLines);

   typedef ColumnIterator::DNASet DNASet;
   typedef ColumnIterator::ColumnMap ColumnMap;
   typedef PositionCache::IntervalSet IntervalSet;
   
   typedef std::pair<const Sequence*, hal_size_t> SeqIndex;
   typedef std::map<SeqIndex, PositionCache*> PositionMap;
   
protected: 
   
   ColumnIteratorConstPtr _colIt;
   std::set<std::string> _missedSet;
   bool _outParalogy;
};

}
#endif
