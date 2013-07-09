/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALWIGGLELIFTOVER_H
#define _HALWIGGLELIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "halWiggleScanner.h"

namespace hal {

class WiggleLiftover : public WiggleScanner
{
public:
   
   WiggleLiftover();
   virtual ~WiggleLiftover();

   void convert(AlignmentConstPtr alignment,
                const Genome* srcGenome,
                std::istream* inputFile,
                const Genome* tgtGenome,
                std::ostream* outputFile,
                bool traverseDupes = true,
                bool unique = false);

   virtual void visitLine();
   virtual void visitHeader();
                      
protected: 

   AlignmentConstPtr _alignment;
   std::istream* _inStream;
   std::ostream* _outStream;
   bool _traverseDupes;
   bool _unique;
      
   const Genome* _srcGenome;
   const Genome* _tgtGenome;
   const Sequence* _srcSequence;
   std::set<const Genome*> _tgtSet;

};

}
#endif
