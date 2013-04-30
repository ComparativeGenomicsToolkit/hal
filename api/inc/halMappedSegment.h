/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAPPEDSEGMENT_H
#define _HALMAPPEDSEGMENT_H

#include "halDefs.h"
#include "halSlicedSegment.h"

namespace hal {

/** 
 * Interface for a mapped segement.  A mapped segment keeps track of a 
 * homologous region in another genome (from which it was mapped).  
 * Mapped segments are used to keep
 * pairwise alignment fragments across the tree as an alternative to the column
 * iterator. 
 * 
 */
class MappedSegment : public virtual SlicedSegment
{
public:

   /** Get the original segment from which this segment was mapped */
   virtual SlicedSegmentConstPtr getSource() const = 0;

protected:
   friend class counted_ptr<MappedSegment>;
   friend class counted_ptr<const MappedSegment>;
   virtual ~MappedSegment() = 0;
};

inline MappedSegment::~MappedSegment() {}

}
#endif
