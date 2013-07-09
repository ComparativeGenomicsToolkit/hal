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
#include "hal.h"

namespace hal {

class WiggleLiftover
{
public:
   
   WiggleLiftover();
   virtual ~WiggleLiftover();

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

   AlignmentConstPtr _alignment;
   std::ostream* _outBedStream;
   bool _traverseDupes;
      
   const Genome* _srcGenome;
   const Genome* _tgtGenome;
   const Sequence* _srcSequence;
   std::set<const Genome*> _tgtSet;

};

}
#endif
