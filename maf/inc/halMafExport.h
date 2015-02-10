/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFEXPORT_H
#define _HALMAFEXPORT_H

#include <iostream>
#include <string>
#include <set>
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
                                 const std::set<const Genome*>& targets);

   // Convert all columns in the leaf genomes to MAF. Each column is
   // reported exactly once regardless of the unique setting, although
   // this may change in the future. Likewise, maxRefGap has no
   // effect, although noDupes will work.
   void convertEntireAlignment(std::ostream& mafStream,
                               AlignmentConstPtr alignment);

   void setMaxRefGap(hal_size_t maxRefGap);
   void setNoDupes(bool noDupes);
   void setNoAncestors(bool noAncestors);
   void setUcscNames(bool ucscNames);
   void setUnique(bool unique);
   void setAppend(bool append);
   void setMaxBlockLength(hal_index_t maxLength);
   void setPrintTree(bool printTree);
   void setOnlyOrthologs(bool onlyOrthologs);

protected:

   void writeHeader();

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _mafStream;
   MafBlock _mafBlock;
   hal_size_t _maxRefGap;
   bool _noDupes;
   bool _noAncestors;
   bool _ucscNames;
   bool _unique;
   bool _append;
   bool _printTree;
   bool _onlyOrthologs;
};

}

#endif
