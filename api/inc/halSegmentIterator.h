/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSlicedSegment.h"
#include "halMappedSegmentContainers.h"

namespace hal {

/** 
 * Interface for general segment iterator.  Common functionality
 * for top and bottom iterators.  A segment iterator implements 
 * the segment as well as the iterator interface. 
 */
class  SegmentIterator : public virtual SlicedSegment
{
public:
    /* constructor */
    SegmentIterator(hal_offset_t startOffset = 0, 
                    hal_offset_t endOffset = 0,
                    bool reversed = false);
    
   /** Destructor */
    virtual ~SegmentIterator() {
   }

    /** move iterator one segment to the left 
    * @param leftCutoff If the left segment contains the 
    * genome position (DNA coordinate) specified by leftCutoff
    * the segment iterator is sliced to begin at this position 
    * (if the iterator is reversed, it moves in opposite direction) */
    virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX);

   /** move iterator one segment to the right 
    * @param rightCutoff If the right segment contains the 
    * genome position (DNA coordinate) specified by rightCutoff
    * the segment iterator is sliced to end at this position 
    * (if the iterator is reversed, it moves in opposite direction) */
    virtual void toRight(hal_index_t rightCutoff = NULL_INDEX);

   /** move iterator to position of dna site in *genome*
    * @param position index of site in genome
    * @param slice if true, the returned iterator is sliced to correspond
    * to just the single base.  if false, the iterator corresponds to the
    * entire segment. 
    * NOTE*** this function requires up to log2(N) time in current hdf5 imp.
    * though it should be faster on average*/
    virtual void toSite(hal_index_t position, bool slice = true);

    /** has the iterator reach the end of the traversal in the direction of
     * movement? */
    bool atEnd() const {
        if (not _reversed) {
            return ((hal_size_t)getArrayIndex() >= getNumSegmentsInGenome());
        } else {
            return (getArrayIndex() < 0);
        }
    }
    
    // FIXME: document or change way getting segment works
    virtual Segment* getSegment() = 0;
   virtual const Segment* getSegment() const = 0;

   // SEGMENT INTERFACE
   virtual void setArrayIndex(Genome* genome, 
                              hal_index_t arrayIndex);
   virtual const Genome* getGenome() const;
   virtual Genome* getGenome();
   virtual const Sequence* getSequence() const;
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
     MappedSegmentSet& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath = NULL,
     bool doDupes = true,
     hal_size_t minLength = 0,
     const Genome *coalescenceLimit = NULL,
     const Genome *mrca = NULL) const;

    virtual void print(std::ostream& os) const;


    
   // SLICED SEGMENT INTERFACE 
   virtual void toReverse();
   virtual void toReverseInPlace();
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset );
   virtual bool getReversed() const;
protected:
   virtual bool inRange() const;
   virtual hal_size_t getNumSegmentsInGenome() const;
   
   hal_offset_t _startOffset;
   hal_offset_t _endOffset;
   bool _reversed;
};

HAL_FORWARD_DEC_MUTABLE_CLASS(SegmentIterator)

inline bool SegmentIterator::inRange() const
{
  return getArrayIndex() >= 0 && getArrayIndex() < 
     (hal_index_t)getNumSegmentsInGenome();
}

inline bool operator<(SegmentIteratorPtr segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt->leftOf(genomePos);
}

inline bool operator>(SegmentIteratorPtr segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt->rightOf(genomePos);
}

}
#endif
