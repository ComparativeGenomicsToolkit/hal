/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDTOPSEGMENTITERATOR_H
#define _HALGAPPEDTOPSEGMENTITERATOR_H

#include <iostream>
#include "halDefs.h"
#include "halTopSegmentIterator.h"
#include "halGappedSegmentIterator.h"

namespace hal {

/** 
 * Interface for Gapped Top Segment iterator.  Only used internally
 * and probably shouldn't be in public interface.
 */
class GappedTopSegmentIterator : virtual public GappedSegmentIterator
{
public:
    /** constructor */
    GappedTopSegmentIterator(TopSegmentIteratorConstPtr left,
                             hal_size_t gapThreshold,
                             bool atomic);

    
   /** Destructor */
    virtual ~GappedTopSegmentIterator() {
    }
    
   /** Return a copy of the iterator */
   virtual GappedTopSegmentIteratorPtr copy();

   /** Return a copy of the iterator */
   virtual GappedTopSegmentIteratorConstPtr copy() const;

   /** Copy another iterator into the current iterator (more efficient
    * than above methods since no new iterators are created */
   virtual void copy(GappedTopSegmentIteratorConstPtr ts) const;

   /** Move to child of given iterator 
    * @param bs Move to child of this iterator */
   virtual void toChild(GappedBottomSegmentIteratorConstPtr bs) const;

   /** Test equality with other iterator 
    * @param other */
   virtual bool equals(GappedTopSegmentIteratorConstPtr other) const;

   /** Test if iterator abuts other iterator */
   virtual bool adjacentTo(GappedTopSegmentIteratorConstPtr other) const;

   /** Test if iterator has a parent */
   virtual bool hasParent() const;

   /** Test if iterator has a duplicate */
   virtual bool hasNextParalogy() const;

   /** Move to next paralogy */
   virtual void toNextParalogy() const;

   /** Test if in reverse orientationt with respect to parent */
   virtual bool getParentReversed() const;

   /** Return the leftmost segment of the iterator
    * (note that moving the returned iterator will corrupt the 
    * current gapped iterator.  this is a bug) */
   virtual TopSegmentIteratorConstPtr getLeft() const;

   /** Return the rightmost segment of the iterator
    * (note that moving the returned iterator will corrupt the 
    * current gapped iterator.  this is a bug) */
   virtual TopSegmentIteratorConstPtr getRight() const;

   /** Reset the gapped iterator.
    * @param ts This will be the left segment of the current iterator. The 
    * right segment will be extended as far as possible */
   virtual void setLeft(TopSegmentIteratorConstPtr ts) const;

   /** For every set of paralogous top segments in a given genome, we identify a
    * unique segment as the canonical reference.  This is the segment that 
    * will be traversed when disabling duplication edges when mapping via 
    * the column iterator or mappedSegment interface.  It is also the segment
    * that is connected from its parent's down edge.*/
   virtual bool isCanonicalParalog() const;

    /* get the current segment */
    virtual SegmentPtr getSegment();

    /* get the current segment */
    virtual SegmentConstPtr getSegment() const;
    
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
     MappedSegmentConstSet& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath,
     bool doDupes,
     hal_size_t minLength,
     const Genome *coalescenceLimit,
     const Genome *mrca) const;
   virtual void print(std::ostream& os) const;

   // SEGMENT ITERATOR INTERFACE
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   virtual void toReverse() const;
   virtual void toReverseInPlace() const;
   virtual void toSite(hal_index_t position, bool slice = true) const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset,
                      hal_offset_t endOffset) const;
   virtual bool getReversed() const;

   // GAPPED SEGMENT ITERATOR INTERFACE
   virtual hal_size_t getGapThreshold() const;
   virtual bool getAtomic() const;
   virtual hal_size_t getChildIndex() const;
   virtual hal_size_t getNumSegments() const;
   virtual hal_size_t getNumGaps() const;
   virtual hal_size_t getNumGapBases() const;
   virtual hal_index_t getLeftArrayIndex() const;
   virtual hal_index_t getRightArrayIndex() const;

private:
   
   bool compatible(TopSegmentIteratorConstPtr left,
                   TopSegmentIteratorConstPtr right) const;

   void extendRight() const;
   void extendLeft() const;

   void toLeftNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toRightNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toLeftNextUngapped(TopSegmentIteratorConstPtr ts) const;
   void toRightNextUngapped(TopSegmentIteratorConstPtr ts) const;
   
   // keep convention of other iterators where const-ness only applies
   // to the database and not the iterator...
   mutable TopSegmentIteratorConstPtr _left;
   mutable TopSegmentIteratorConstPtr _right;
   mutable BottomSegmentIteratorConstPtr _leftParent;
   mutable BottomSegmentIteratorConstPtr _rightParent;
   mutable TopSegmentIteratorConstPtr _leftDup;
   mutable TopSegmentIteratorConstPtr _rightDup;
   mutable TopSegmentIteratorConstPtr _temp;
   mutable TopSegmentIteratorConstPtr _temp2;
   mutable hal_size_t _childIndex;
   mutable hal_size_t _gapThreshold;
   mutable bool _atomic;
  
};


inline bool operator==(GappedTopSegmentIteratorConstPtr p1,
                       GappedTopSegmentIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(GappedTopSegmentIteratorConstPtr p1,
                       GappedTopSegmentIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}

#endif
