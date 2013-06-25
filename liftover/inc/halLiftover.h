/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLIFTOVER_H
#define _HALLIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "halBedScanner.h"

namespace hal {

class Liftover : public BedScanner
{
public:
   
   Liftover();
   virtual ~Liftover();

   void convert(AlignmentConstPtr alignment,
                const Genome* srcGenome,
                std::istream* inputFile,
                const Genome* tgtGenome,
                std::ostream* outputFile,
                int inBedVersion = -1,
                int outBedVersion = -1,
                bool addExtraColumns = false,
                bool traverseDupes = true,
                bool outPSL = false);
                   
protected:

   typedef std::list<BedLine> BedList;

   virtual void visitBegin();
   virtual void visitLine();
   virtual void visitEOF();
   virtual void writeLineResults();
   virtual void assignBlocksToIntervals();
   virtual void writeBlocksAsIntervals();
   virtual void cleanResults();
   virtual void liftBlockIntervals();
   virtual void mergeIntervals();
   virtual void liftInterval(BedList& mappedBedLines) = 0;
   
protected: 

   AlignmentConstPtr _alignment;
   std::ostream* _outBedStream;
   bool _addExtraColumns;
   bool _traverseDupes;
   BedList _outBedLines;
   int _inBedVersion;
   int _outBedVersion;
   bool _outPSL;
   
   BedList _mappedBlocks;
   
   const Genome* _srcGenome;
   const Genome* _tgtGenome;
   const Sequence* _srcSequence;
   std::set<const Genome*> _tgtSet;

   ColumnIteratorConstPtr _colIt;
   std::set<std::string> _missedSet;
};

}
#endif
