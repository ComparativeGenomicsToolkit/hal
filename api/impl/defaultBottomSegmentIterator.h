/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTBOTTOMSEGMENTITERATOR_H
#define _DEFAULTBOTTOMSEGMENTITERATOR_H

#include "halBottomSegmentIterator.h"
#include "defaultSegmentIterator.h"

namespace hal {

class DefaultBottomSegmentIterator : public DefaultSegmentIterator,
                                     virtual public BottomSegmentIterator
{
public:
   
   DefaultBottomSegmentIterator(BottomSegment* segment,
                                hal_size_t startOffset = 0, 
                                hal_size_t endOffset = 0,
                                bool inverted = false);
   virtual ~DefaultBottomSegmentIterator();
   
   // SEGMENT INTERFACE OVERRIDE
   virtual void print(std::ostream& os) const;

   // BOTTOM SEGMENT INTERFACE
   virtual hal_size_t getNumChildren() const;
   virtual hal_index_t getChildIndex(hal_size_t i) const;
   virtual hal_index_t getChildIndexG(const Genome* childGenome) const;
   virtual bool hasChild(hal_size_t child) const;
   virtual bool hasChildG(const Genome* childGenome) const;
   virtual void setChildIndex(hal_size_t i, hal_index_t childIndex);
   virtual bool getChildReversed(hal_size_t i) const;
   virtual void setChildReversed(hal_size_t child, bool isReversed);
   virtual hal_index_t getTopParseIndex() const;
   virtual void setTopParseIndex(hal_index_t parseIndex);
   virtual hal_offset_t getTopParseOffset() const;
   virtual bool hasParseUp() const;
   virtual hal_index_t getLeftChildIndex(hal_size_t i) const;
   virtual hal_index_t getRightChildIndex(hal_size_t i) const;

   // BOTTOM SEGMENT ITERATOR INTERFACE
   virtual BottomSegmentIteratorPtr copy();
   virtual BottomSegmentIteratorConstPtr copy() const;
   virtual void copy(BottomSegmentIteratorConstPtr bs) const;
   virtual void toParent(TopSegmentIteratorConstPtr ts) const; 
   virtual void toParseDown(TopSegmentIteratorConstPtr ts) const;
   virtual BottomSegment* getBottomSegment();
   virtual const BottomSegment* getBottomSegment() const;
   virtual bool equals(BottomSegmentIteratorConstPtr other) const;

protected:
   virtual SegmentPtr getSegment();
   virtual SegmentConstPtr getSegment() const;
   virtual hal_size_t getNumSegmentsInGenome() const;

protected:

   BottomSegmentPtr _bottomSegment;

};

}
#endif
