/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTTOPSEGMENTITERATOR_H
#define _DEFAULTTOPSEGMENTITERATOR_H

#include "halTopSegmentIterator.h"
#include "defaultSegmentIterator.h"

namespace hal {


class DefaultTopSegmentIterator : virtual public DefaultSegmentIterator, 
                                  virtual public TopSegmentIterator
{
public:
   DefaultTopSegmentIterator(TopSegment* topSegment,
                             hal_offset_t startOffset = 0, 
                             hal_offset_t endOffset = 0,
                             bool inverted = false);
   virtual ~DefaultTopSegmentIterator();
   
   // SEGMENT INTERFACE OVERRIDE
   virtual void print(std::ostream& os) const;

   // TOP SEGMENT INTERFACE
   virtual hal_index_t getParentIndex() const;
   virtual bool hasParent() const;
   virtual void setParentIndex(hal_index_t parIdx);
   virtual bool getParentReversed() const;
   virtual void setParentReversed(bool isReversed);
   virtual hal_index_t getBottomParseIndex() const;
   virtual void setBottomParseIndex(hal_index_t botParseIdx);
   virtual hal_offset_t getBottomParseOffset() const;
   virtual bool hasParseDown() const;
   virtual hal_index_t getNextParalogyIndex() const;
   virtual bool hasNextParalogy() const;
   virtual void setNextParalogyIndex(hal_index_t parIdx);
   virtual hal_index_t getLeftParentIndex() const;
   virtual hal_index_t getRightParentIndex() const;
   virtual bool isCanonicalParalog() const;


   // TOP SEGMENT ITERATOR INTERFACE 
   virtual TopSegmentIteratorPtr copy();
   virtual TopSegmentIteratorConstPtr copy() const;
   virtual void copy(TopSegmentIteratorConstPtr ts) const;
   virtual void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const;
   virtual void toChildG(BottomSegmentIteratorConstPtr bs, 
                         const Genome* childGenome) const;
   virtual void toParseUp(BottomSegmentIteratorConstPtr bs) const;
   virtual TopSegment* getTopSegment();
   virtual const TopSegment* getTopSegment() const;
   virtual bool equals(TopSegmentIteratorConstPtr other) const;
   virtual void toNextParalogy() const;

protected:
   virtual SegmentPtr getSegment();
   virtual SegmentConstPtr getSegment() const;
   virtual hal_size_t getNumSegmentsInGenome() const;

protected:

   TopSegmentPtr _topSegment;
};




}
#endif
