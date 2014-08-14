/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTMAPPEDSEGMENT_H
#define _DEFAULTMAPPEDSEGMENT_H

#include <list>
#include "halMappedSegment.h"
#include "defaultSegmentIterator.h"

namespace hal {

HAL_FORWARD_DEC_CLASS(DefaultMappedSegment)

// note it would be nice to extend DefaultSegmentIterator but the 
// crappy smart pointer interface makes it impossible to use "this"
// as a parameter to lots of api functions.  simplest for now just
// to contain a pair of segment iterators and wrap up all the interface
// methods. 
class DefaultMappedSegment : virtual public MappedSegment
{
public:
   
   virtual ~DefaultMappedSegment();
   
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
     std::set<MappedSegmentConstPtr>& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath,
     bool doDupes,
     hal_size_t minLength,
     const Genome *coalescenceLimit,
     const Genome *mrca) const;
   virtual void print(std::ostream& os) const;

   // SLICED SEGMENT INTERFACE 
   virtual void toReverse() const;
   virtual void toReverseInPlace() const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset ) const;
   virtual bool getReversed() const;

   // MAPPED SEGMENT INTERFACE 
   virtual SlicedSegmentConstPtr getSource() const;
   virtual bool lessThan(const MappedSegmentConstPtr& other) const;
   virtual bool lessThanBySource(const MappedSegmentConstPtr& other) const;
   virtual bool equals(const MappedSegmentConstPtr& other) const;
   virtual void flip() const;
   virtual void fullReverse() const;
   virtual MappedSegmentConstPtr copy() const;
   virtual bool canMergeRightWith(
     const MappedSegmentConstPtr& next,
     const std::set<hal_index_t>* cutSet,
     const std::set<hal_index_t>* sourceCutSet) const;


   // INTERNAL METHODS
   static hal_size_t map(const DefaultSegmentIterator* source,
                         std::set<MappedSegmentConstPtr>& results,
                         const Genome* tgtGenome,
                         const std::set<const Genome*>* genomesOnPath,
                         bool doDupes,
                         hal_size_t minLength,
                         const Genome *coalescenceLimit,
                         const Genome *mrca);


protected:
   friend class counted_ptr<DefaultMappedSegment>;
   friend class counted_ptr<const DefaultMappedSegment>;

protected:

   DefaultMappedSegment(SegmentIteratorConstPtr source,
                        SegmentIteratorConstPtr target);

   static 
   int fastComp(const DefaultSegmentIteratorConstPtr& s1, 
                const DefaultSegmentIteratorConstPtr& s2);

   static 
   int boundComp(const DefaultSegmentIteratorConstPtr& s1, 
                const DefaultSegmentIteratorConstPtr& s2);

   static 
   int slowComp(const DefaultSegmentIteratorConstPtr& s1, 
                const DefaultSegmentIteratorConstPtr& s2);
   
   enum OverlapCat { Same, Disjoint, AContainsB, BContainsA,
                     AOverlapsLeftOfB, BOverlapsLeftOfA };
 
   static
   OverlapCat slowOverlap(const SlicedSegmentConstPtr& s1, 
                          const SlicedSegmentConstPtr& s2);

   static
   void getOverlapBounds(const DefaultMappedSegmentConstPtr& seg, 
                         std::set<MappedSegmentConstPtr>& results, 
                         std::set<MappedSegmentConstPtr>::iterator& leftBound, 
                         std::set<MappedSegmentConstPtr>::iterator& rightBound);

   static 
   void clipAagainstB(MappedSegmentConstPtr segA,
                      MappedSegmentConstPtr segB,
                      OverlapCat overlapCat,
                      std::vector<MappedSegmentConstPtr>& clippedSegs);
   static 
   void insertAndBreakOverlaps(DefaultMappedSegmentConstPtr seg,
                               std::set<MappedSegmentConstPtr>& results);
   
   // Map a segment to all segments that share any homology in or below
   // the given "coalescence limit" genome (not just those that share
   // homology in the MRCA of the source and target genomes).
   static hal_size_t mapIncludingExtraParalogs(
     const Genome* srcGenome,
     std::list<DefaultMappedSegmentConstPtr>& input,
     std::list<DefaultMappedSegmentConstPtr>& results,
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
     std::list<DefaultMappedSegmentConstPtr>& input,
     std::list<DefaultMappedSegmentConstPtr>& results,
     const std::set<std::string>& namesOnPath,
     const Genome* coalescenceLimit,
     hal_size_t minLength);

   // Map the input segments up until reaching the target genome. If the
   // target genome is below the source genome, fail miserably.
   // Destructive to any data in the input or results list.
   static hal_size_t mapRecursiveUp(
     std::list<DefaultMappedSegmentConstPtr>& input,
     std::list<DefaultMappedSegmentConstPtr>& results,
     const Genome* tgtGenome,
     hal_size_t minLength);

   // Map the input segments down until reaching the target genome. If the
   // target genome is above the source genome, fail miserably.
   // Destructive to any data in the input or results list.
   static hal_size_t mapRecursiveDown(
     std::list<DefaultMappedSegmentConstPtr>& input,
     std::list<DefaultMappedSegmentConstPtr>& results,
     const Genome* tgtGenome,
     const std::set<std::string>& namesOnPath,
     bool doDupes,
     hal_size_t minLength);

   static 
   hal_size_t mapRecursive(const Genome* prevGenome,
                           std::list<DefaultMappedSegmentConstPtr>& input,
                           std::list<DefaultMappedSegmentConstPtr>& results,
                           const Genome* tgtGenome,
                           const std::set<std::string>& namesOnPath,
                           bool doDupes,
                           hal_size_t minLength);
   static 
   hal_size_t mapUp(DefaultMappedSegmentConstPtr mappedSeg, 
                    std::list<DefaultMappedSegmentConstPtr>& results,
                    bool doDupes,
                    hal_size_t minLength);
   static 
   hal_size_t mapDown(DefaultMappedSegmentConstPtr mappedSeg, 
                      hal_size_t childIndex,
                      std::list<DefaultMappedSegmentConstPtr>& results,
                      hal_size_t minLength);
   static 
   hal_size_t mapSelf(DefaultMappedSegmentConstPtr mappedSeg, 
                      std::list<DefaultMappedSegmentConstPtr>& results,
                    hal_size_t minLength);
   
   TopSegmentIteratorConstPtr targetAsTop() const;
   BottomSegmentIteratorConstPtr targetAsBottom() const;
   TopSegmentIteratorConstPtr sourceAsTop() const;
   BottomSegmentIteratorConstPtr sourceAsBottom() const;
   SegmentIteratorConstPtr sourceCopy() const;

  struct LessSource {
      bool operator()(const DefaultMappedSegmentConstPtr& ms1,
                      const DefaultMappedSegmentConstPtr& ms2) const;
        };
   
   struct Less {
      bool operator()(const DefaultMappedSegmentConstPtr& ms1,
                      const DefaultMappedSegmentConstPtr& ms2) const;
        };
   
   struct EqualTo {
      bool operator()(const DefaultMappedSegmentConstPtr& ms1,
                      const DefaultMappedSegmentConstPtr& ms2) const;
        };
   
protected:

   mutable DefaultSegmentIteratorConstPtr _source;
   mutable DefaultSegmentIteratorConstPtr _target;
};

inline TopSegmentIteratorConstPtr DefaultMappedSegment::targetAsTop() const
{
  return _target.downCast<TopSegmentIteratorConstPtr>();
}

inline 
BottomSegmentIteratorConstPtr DefaultMappedSegment::targetAsBottom() const
{
  return _target.downCast<BottomSegmentIteratorConstPtr>();
}

inline TopSegmentIteratorConstPtr DefaultMappedSegment::sourceAsTop() const
{
  return _source.downCast<TopSegmentIteratorConstPtr>();
}

inline 
BottomSegmentIteratorConstPtr DefaultMappedSegment::sourceAsBottom() const
{
  return _source.downCast<BottomSegmentIteratorConstPtr>();
}

inline SegmentIteratorConstPtr DefaultMappedSegment::sourceCopy() const
{
  if (_source->isTop())
  {
    return sourceAsTop()->copy();
  }
  return sourceAsBottom()->copy();
}

inline bool DefaultMappedSegment::LessSource::operator()(
  const DefaultMappedSegmentConstPtr& ms1,
  const DefaultMappedSegmentConstPtr& ms2) const
{
  return ms1->lessThanBySource(ms2); 
};

inline bool DefaultMappedSegment::Less::operator()(
  const DefaultMappedSegmentConstPtr& ms1,
  const DefaultMappedSegmentConstPtr& ms2) const
{
  return ms1->lessThan(ms2); 
};

inline bool DefaultMappedSegment::EqualTo::operator()(
  const DefaultMappedSegmentConstPtr& ms1,
  const DefaultMappedSegmentConstPtr& ms2) const
{
  return ms1->equals(ms2); 
};



}
#endif
