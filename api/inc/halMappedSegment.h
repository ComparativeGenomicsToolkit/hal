/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAPPEDSEGMENT_H
#define _HALMAPPEDSEGMENT_H

#include <list>
#include "halDefs.h"
#include "halSlicedSegment.h"
#include "halSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include "halBottomSegmentIterator.h"

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

   /** Destructor */
    virtual ~MappedSegment() {
    }

    /** Get the original segment from which this segment was mapped */
    virtual SlicedSegment* getSource() const {
        return _source.get();
    }

    /** Get the original segment from which this segment was mapped */
    virtual SlicedSegmentPtr getSourcePtr() const {
        return _source;
    }

   /** Comparison used to store in stl sets and maps.  We sort based
    * on the coordinate of the mapped segemnt's (target) interval as the primary
    * index and the target genome as the secondary index.  */
   virtual bool lessThan(const MappedSegment* other) const;

   /** Comparison used to store in STL containers.  We sort based
    * on the coordinate of the *Source* genome interval as the primary
    * index and the target genome as the secondary index.  */
   virtual bool lessThanBySource(const MappedSegment* other) const;

   /** Comparison used to determine uniqueness in lists.  Tests lessThan
    * in both directions.  Note that equality is the same regardless of 
    * whether or not we use the source segment as our primary index. */
   virtual bool equals(const MappedSegment* other) const;

   /** Flip the mapping direction.  This segment becomes the source, and
    * the source becomes this.*/
   virtual void flip();

   /** Reverse both segments.  Also swap their start and end offsets. 
    * Note that toReverse() does not reverse the source segment.*/
   virtual void fullReverse();
   
   /** Return of a copy of the mapped segment */
   virtual MappedSegment* clone() const;

   /** Test if mapped segment can be merged to the right with input 
    * segment.  will return false if the right coordinate of this is in
    * either (optional) cutSet.*/
   virtual bool canMergeRightWith(
     const MappedSegmentPtr& next,
     const std::set<hal_index_t>* cutSet = NULL,
     const std::set<hal_index_t>* sourceCutSet = NULL) const;

    // FIXME: are non-ptr functors needed?
   /** Functor for sorted STL containers, sorting by origin as primary 
    * index */
   struct LessSource { bool operator()(const MappedSegment& ms1,
                                       const MappedSegment& ms2) const {
     return ms1.lessThanBySource(&ms2); }
   };

   /** Functor for sorted STL containers, sorting by origin as primary 
    * index */
   struct LessSourcePtr { bool operator()(const MappedSegmentPtr& ms1,
                                          const MappedSegmentPtr& ms2) const {
       return ms1->lessThanBySource(ms2.get());
   }
   };

   /** Functor for sorted STL containers, sorting  by target as primary 
    * index */
   struct Less { bool operator()(const MappedSegment& ms1,
                                 const MappedSegment& ms2) const {
     return ms1.lessThan(&ms2); }
   };
 
   /** Functor for sorted STL containers, sorting  by target as primary 
    * index */
   struct LessPtr { bool operator()(const MappedSegmentPtr& ms1,
                                 const MappedSegmentPtr& ms2) const {
       return ms1->lessThan(ms2.get()); }
   };
 
   /** Functor for STL sorted lists to test for uniqueness */
   struct EqualTo { bool operator()(const MappedSegment& ms1,
                                    const MappedSegment& ms2) const {
     return ms1.equals(&ms2); }
   };

   /** Functor for STL sorted lists to test for uniqueness */
   struct EqualToPtr { bool operator()(const MappedSegmentPtr& ms1,
                                    const MappedSegmentPtr& ms2) const {
       return ms1->equals(ms2.get()); }
   };

   // NEEDS TO BE ADDED TO SEGMENT INTERFACE
    virtual void print(std::ostream& os) const;

   // SEGMENT INTERFACE
   virtual void setArrayIndex( Genome* genome, 
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

   // SLICED SEGMENT INTERFACE 
   virtual void toReverse();
   virtual void toReverseInPlace();
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset,
                      hal_offset_t endOffset );
   virtual bool getReversed() const;

   // INTERNAL METHODS
   static hal_size_t map(const SegmentIterator* source,
                         MappedSegmentSet& results,
                         const Genome* tgtGenome,
                         const std::set<const Genome*>* genomesOnPath,
                         bool doDupes,
                         hal_size_t minLength,
                         const Genome *coalescenceLimit,
                         const Genome *mrca);
private:

   MappedSegment(SegmentIteratorPtr source,
                 SegmentIteratorPtr target);

   static 
   int fastComp(const SegmentIteratorPtr& s1, 
                const SegmentIteratorPtr& s2);

   static 
   int boundComp(const SegmentIteratorPtr& s1, 
                const SegmentIteratorPtr& s2);

   static 
   int slowComp(const SegmentIteratorPtr& s1, 
                const SegmentIteratorPtr& s2);
   
   enum OverlapCat { Same, Disjoint, AContainsB, BContainsA,
                     AOverlapsLeftOfB, BOverlapsLeftOfA };
 
   static
   OverlapCat slowOverlap(const SlicedSegment* s1, 
                          const SlicedSegment* s2);

   static
   void getOverlapBounds(MappedSegmentPtr& seg, 
                         MappedSegmentSet& results, 
                         MappedSegmentSet::iterator& leftBound, 
                         MappedSegmentSet::iterator& rightBound);

   static 
   void clipAagainstB(MappedSegmentPtr segA,
                      MappedSegmentPtr segB,
                      OverlapCat overlapCat,
                      std::vector<MappedSegmentPtr>& clippedSegs);
   static 
   void insertAndBreakOverlaps(MappedSegmentPtr seg,
                               MappedSegmentSet& results);
   
   // Map a segment to all segments that share any homology in or below
   // the given "coalescence limit" genome (not just those that share
   // homology in the MRCA of the source and target genomes).
   static hal_size_t mapIncludingExtraParalogs(
     const Genome* srcGenome,
     std::list<MappedSegmentPtr>& input,
     std::list<MappedSegmentPtr>& results,
     const std::set<std::string>& namesOnPath,
     const Genome* tgtGenome,
     const Genome* mrca,
     const Genome *coalescenceLimit,
     bool doDupes,
     hal_size_t minLength);

   // Map all segments from the input to any segments in the same genome
   // that coalesce in or before the given "coalescence limit" genome.
   // Destructive to any data in the input list.
   static hal_size_t mapRecursiveParalogies(
     const Genome *srcGenome,
     std::list<MappedSegmentPtr>& input,
     std::list<MappedSegmentPtr>& results,
     const std::set<std::string>& namesOnPath,
     const Genome* coalescenceLimit,
     hal_size_t minLength);

   // Map the input segments up until reaching the target genome. If the
   // target genome is below the source genome, fail miserably.
   // Destructive to any data in the input or results list.
   static hal_size_t mapRecursiveUp(
     std::list<MappedSegmentPtr>& input,
     std::list<MappedSegmentPtr>& results,
     const Genome* tgtGenome,
     hal_size_t minLength);

   // Map the input segments down until reaching the target genome. If the
   // target genome is above the source genome, fail miserably.
   // Destructive to any data in the input or results list.
   static hal_size_t mapRecursiveDown(
     std::list<MappedSegmentPtr>& input,
     std::list<MappedSegmentPtr>& results,
     const Genome* tgtGenome,
     const std::set<std::string>& namesOnPath,
     bool doDupes,
     hal_size_t minLength);

   static 
   hal_size_t mapRecursive(const Genome* prevGenome,
                           std::list<MappedSegmentPtr>& input,
                           std::list<MappedSegmentPtr>& results,
                           const Genome* tgtGenome,
                           const std::set<std::string>& namesOnPath,
                           bool doDupes,
                           hal_size_t minLength);
   static 
   hal_size_t mapUp(MappedSegmentPtr mappedSeg, 
                    std::list<MappedSegmentPtr>& results,
                    bool doDupes,
                    hal_size_t minLength);
   static 
   hal_size_t mapDown(MappedSegmentPtr mappedSeg, 
                      hal_size_t childIndex,
                      std::list<MappedSegmentPtr>& results,
                      hal_size_t minLength);
   static 
   hal_size_t mapSelf(MappedSegmentPtr mappedSeg, 
                      std::list<MappedSegmentPtr>& results,
                      hal_size_t minLength);
   
   TopSegmentIteratorPtr targetAsTop() const;
   BottomSegmentIteratorPtr targetAsBottom() const;
   TopSegmentIteratorPtr sourceAsTop() const;
   BottomSegmentIteratorPtr sourceAsBottom() const;
   SegmentIteratorPtr sourceclone() const;

   
private:

   SegmentIteratorPtr _source;
   SegmentIteratorPtr _target;
};
}

inline hal::TopSegmentIteratorPtr hal::MappedSegment::targetAsTop() const
{
    return std::dynamic_pointer_cast<TopSegmentIterator>(_target);
}

inline 
hal::BottomSegmentIteratorPtr hal::MappedSegment::targetAsBottom() const
{
    return std::dynamic_pointer_cast<BottomSegmentIterator>(_target);
}

inline hal::TopSegmentIteratorPtr hal::MappedSegment::sourceAsTop() const
{
    return std::dynamic_pointer_cast<TopSegmentIterator>(_source);
}

inline 
hal::BottomSegmentIteratorPtr hal::MappedSegment::sourceAsBottom() const
{
    return std::dynamic_pointer_cast<BottomSegmentIterator>(_source);
}

inline hal::SegmentIteratorPtr hal::MappedSegment::sourceclone() const
{
  if (_source->isTop()) {
    return sourceAsTop()->clone();
  } else {
      return sourceAsBottom()->clone();
  }
}

#endif
