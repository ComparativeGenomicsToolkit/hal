/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBRANCHMUTATIONS_H
#define _HALBRANCHMUTATIONS_H

#include <iostream>
#include <string>
#include <map>
#include "hal.h"

namespace hal {

class BranchMutations
{
public:

   BranchMutations();
   virtual ~BranchMutations();

   void printCsv(std::ostream& outStream) const;
   void analyzeBranch(AlignmentConstPtr alignment,
                      hal_size_t gapThreshold,
                      std::ostream* refBedStream,
                      std::ostream* delBedStream,
                      std::ostream* snpBedStream,
                      const Genome* reference,
                      hal_index_t startPosition,
                      hal_size_t length);

protected:

   void writeRearrangement();
   void writeSubstitutions();
   void writeGapInsertions();
   void writeGapDeletion();

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _refStream;
   std::ostream* _delStream;
   std::ostream* _snpStream;
   hal_size_t _maxGap;
   const SegmentedSequence* _reference;
   hal_size_t _start;
   hal_size_t _length;

   RearrangementPtr _rearrangement;
};

}

#endif
