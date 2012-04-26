/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for general segment iterator.  Common functionality
 * for top and bottom iterators.  A segment iterator implements 
 * the segment as well as the iterator interface. 
 */
class SegmentIterator 
{
public:

   /** move iterator one segment to the left */
   virtual void toLeft() const = 0;

   /** move iterat one segment to the right */
   virtual void toRight() const = 0;

   /** switch to segment's reverse complement */
   virtual void toReverse() const = 0;

   /** move to the next paralgous segment */
   virtual void toNextParalogy() const = 0;
   
   /** Get the iterator's start offset.  This is used when moving
    * vertically and following the parse index.  Any part of the
    * segment before the start offset is ignored by the iterator */  
   virtual hal_offset_t getStartOffset() const = 0;

   /** Get the iterator's end offset.  This is used when moving
    * vertically and following the parse index.  Any part of the
    * segment after the end offset is ignored by the iterator */  
   virtual hal_offset_t getEndOffset() const = 0;

   /** Get the length of the iterator's segment.  The parse offset's
    * are taken into account.  So this will return the segment's 
    * length minues the sum of the two offsets */
   virtual hal_size_t getLength() const = 0;

   /** Check whether iterator is on segment's reverse complement */
   virtual hal_bool_t getReversed() const = 0;
   
   /** Get the DNA string corresponding to the iterator from the genome 
    * @param outString string into which the results are copied */
   virtual void getString(std::string& outString) const = 0;

protected:
   friend class counted_ptr<SegmentIterator>;
   friend class counted_ptr<const SegmentIterator>;
   virtual ~SegmentIterator() = 0;
};

inline SegmentIterator::~SegmentIterator() {}

}
#endif
