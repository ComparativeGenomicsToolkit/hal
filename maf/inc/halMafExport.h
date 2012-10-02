/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFEXPORT_H
#define _HALMAFEXPORT_H

#include <iostream>
#include <string>
#include <vector>
#include "halMafBlock.h"

namespace hal {

class MafExport
{
public:

   MafExport();

   virtual ~MafExport();

   void convertGenome(std::ostream& mafStream,
                      AlignmentConstPtr alignment,
                      const Genome* genome,
                      hal_index_t startPosition = 0,
                      hal_size_t length = 0,
                      const Genome* root = NULL);

   void convertSequence(std::ostream& mafStream,
                        AlignmentConstPtr alignment,
                        const Sequence* sequence,
                        hal_index_t startPosition = 0,
                        hal_size_t length = 0,
                        const Genome* root = NULL);

   void setMaxRefGap(hal_size_t maxRefGap);
   void setNoDupes(bool noDupes);

protected:

   void writeHeader();

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _mafStream;
   MafBlock _mafBlock;
   hal_size_t _maxRefGap;
   bool _noDupes;
};

}

#endif
