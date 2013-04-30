/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTSLICEDSEGMENT_H
#define _DEFAULTSLICEDSEGMENT_H

#include "halSlicedSegment.h"

namespace hal {


class DefaultSlicedSegment : virtual public SlicedSegment
{
public:
   DefaultSlicedSegment(hal_offset_t startOffset = 0, 
                        hal_offset_t endOffset = 0,
                        bool inverted = false);
   virtual ~DefaultSlicedSegment();
   
   // SEGMENT INTERFACE
   virtual void setArrayIndex(Genome* genome, 
                              hal_index_t arrayIndex);
   virtual void setArrayIndex(const Genome* genome, 
                              hal_index_t arrayIndex) const;
   virtual const Genome* getGenome() const;
   virtual Genome* getGenome();
   virtual const Sequence* getSequence() const;
   virtual Sequence* getSequence();
   virtual hal_index_t getStartPosition() const;
   virtual hal_index_t getEndPosition() const;
   virtual hal_size_t getLength() const;
   virtual void getString(std::string& outString) const;
   virtual void setCoordinates(hal_index_t startPos, hal_size_t length);
   virtual hal_index_t getArrayIndex() const;
   virtual bool leftOf(hal_index_t genomePos) const;
   virtual bool rightOf(hal_index_t genomePos) const;
   virtual bool overlaps(hal_index_t genomePos) const;
   virtual bool isFirst() const;
   virtual bool isLast() const;
   virtual bool isMissingData(double nThreshold) const;
   virtual bool isTop() const;
   virtual hal_size_t getMappedSegments(
     const Genome* tgtGenome,
     std::vector<MappedSegmentConstPtr>& outSegments,
     bool doDupes) const;

   // SLICED SEGMENT INTERFACE 
   virtual void toReverse() const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset ) const;
   virtual bool getReversed() const;

protected:
   virtual SegmentPtr getSegment() = 0;
   virtual SegmentConstPtr getSegment() const = 0;
   virtual bool inRange() const;
   virtual hal_size_t getNumSegmentsInGenome() const = 0;
   
protected:
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable bool _reversed;
};

inline bool DefaultSlicedSegment::inRange() const
{
  return getArrayIndex() >= 0 && getArrayIndex() < 
     (hal_index_t)getNumSegmentsInGenome();
}


}
#endif
