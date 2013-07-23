/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFBED_H
#define _HALMAFBED_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include "halBedScanner.h"
#include "halMafExport.h"

namespace hal {

/** Use the halBedScanner to parse a bed file, running mafExport on each
 * line */
class MafBed : public BedScanner
{
public:

   MafBed(std::ostream& mafStream, AlignmentConstPtr alignment,
          const Genome* refGenome, const Sequence* refSequence,
          hal_index_t refStart, hal_size_t refLength,
          std::set<const Genome*>& targetSet,
          MafExport& mafExport);
   virtual ~MafBed();

   void run(std::istream* bedStream, int bedVersion = -1);

protected: 

   virtual void visitLine();

protected:

   std::ostream& _mafStream;
   AlignmentConstPtr _alignment;
   const Genome* _refGenome;
   const Sequence* _refSequence;
   hal_index_t _refStart;
   hal_size_t _refLength;
   std::set<const Genome*>& _targetSet;
   MafExport& _mafExport;
};

}

#endif
