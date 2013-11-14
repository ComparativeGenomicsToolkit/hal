/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALEXTRACT4D_H
#define _HALEXTRACT4D_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <deque>
#include "halBedScanner.h"

namespace hal {

/** Use the halBedScanner to parse a bed file, extracting 4d codons from
 * each position */
class Extract4d : public BedScanner
{
public:

   Extract4d();
   virtual ~Extract4d();

   void run(const Genome* refGenome,
            std::istream* inBedStream, std::ostream* outBedStream,
            int bedVersion = -1, bool conserved = false);

   static const char CodonPrefixTable[2][8];

protected: 

   virtual void visitLine();

   void extractBed4d();
   void extractConservedBed4d();
   void extractBlocks4d();
   void extractConservedBlocks4d();
   void write();
   

protected:

   std::ostream* _outBedStream;
   const Genome* _refGenome;
   const Sequence* _refSequence;
   std::deque<BedLine> _outBedLines;
   int _bedVersion;
   bool _conserved;
};

}

#endif
