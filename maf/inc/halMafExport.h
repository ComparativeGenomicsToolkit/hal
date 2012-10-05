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

   void convertSegmentedSequence(std::ostream& mafStream,
                                 AlignmentConstPtr alignment,
                                 const SegmentedSequence* seq,
                                 hal_index_t startPosition,
                                 hal_size_t length,
                                 const Genome* root);

   void setMaxRefGap(hal_size_t maxRefGap);
   void setNoDupes(bool noDupes);
   void setNoAncestors(bool noAncestors);


protected:

   void writeHeader();

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _mafStream;
   MafBlock _mafBlock;
   hal_size_t _maxRefGap;
   bool _noDupes;
   bool _noAncestors;
};

}

#endif
