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
                      double nThreshold,
                      std::ostream* refBedStream,
                      std::ostream* parentBedStream,
                      std::ostream* snpBedStream,
                      std::ostream* delBreakBedStream,
                      const Genome* reference,
                      hal_index_t startPosition,
                      hal_size_t length);

   static const std::string inversionBedTag;
   static const std::string insertionBedTag;
   static const std::string deletionBedTag;
   static const std::string deletionBreakBedTag;
   static const std::string transpositionBedTag;
   static const std::string duplicationBedTag;
   static const std::string gapInsertionBedTag;
   static const std::string gapDeletionBedTag;
   static const std::string gapDeletionBreakBedTag;
   static std::string substitutionBedTag(char parent, char child);

protected:

   void writeInsertionOrInversion();
   void writeSubstitutions(TopSegmentIteratorConstPtr first,
                           TopSegmentIteratorConstPtr lastPlusOne);
   void writeGapInsertions();
   void writeDeletion();
   void writeDeletionBreakPoint();
   void writeDuplication();
   void writeHeaders();

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _refStream;
   std::ostream* _parentStream;
   std::ostream* _snpStream;
   std::ostream* _delBreakStream;
   hal_size_t _maxGap;
   double _nThreshold;
   const Genome* _reference;
   const Sequence* _sequence;
   hal_size_t _start;
   hal_size_t _length;
   std::string _refName;
   std::string _parName;

   RearrangementPtr _rearrangement;
   TopSegmentIteratorConstPtr _top;
   BottomSegmentIteratorConstPtr _bottom1, _bottom2;
};

}

#endif
