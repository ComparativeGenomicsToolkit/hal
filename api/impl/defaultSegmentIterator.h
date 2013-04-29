/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTSEGMENTITERATOR_H
#define _DEFAULTSEGMENTITERATOR_H

#include "halSegmentIterator.h"
#include "defaultSlicedSegment.h"

namespace hal {


class DefaultSegmentIterator : public DefaultSlicedSegment,
                               virtual public SegmentIterator
{
public:
   DefaultSegmentIterator(hal_offset_t startOffset = 0, 
                          hal_offset_t endOffset = 0,
                          bool inverted = false);
   virtual ~DefaultSegmentIterator();

   // SEGMENT ITERATOR INTERFACE 
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   virtual void toSite(hal_index_t position, bool slice = true) const;
};

}
#endif
