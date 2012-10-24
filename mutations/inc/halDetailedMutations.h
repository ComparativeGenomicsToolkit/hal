/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDETAILEDMUTATIONS_H
#define _HALDETAILEDMUTATIONS_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"
#include "halMutationsStats.h"

namespace hal {

class DetailedMutations
{
public:

   DetailedMutations();
   virtual ~DetailedMutations();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(hal::AlignmentConstPtr alignment,
                         hal_size_t gapThreshold,
                         std::ostream* snpStream,
                         std::ostream* svStream,
                         const SegmentedSequence* reference,
                         hal_index_t startPosition,
                         hal_size_t length,
                         const std::set<const Genome*>* targets);

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _snpStream;
   std::ostream* _svStream;
   hal_size_t _maxGap;
   ColumnIteratorConstPtr _colIt;

};

}

#endif
