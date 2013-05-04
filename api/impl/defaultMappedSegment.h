/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTMAPPEDSEGMENT_H
#define _DEFAULTMAPPEDSEGMENT_H

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
     std::vector<MappedSegmentConstPtr>& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath,
     bool doDupes) const;

   // SLICED SEGMENT INTERFACE 
   virtual void toReverse() const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset ) const;
   virtual bool getReversed() const;

   // MAPPED SEGMENT INTERFACE 
   virtual SlicedSegmentConstPtr getSource() const;

   // INTERNAL METHODS
   static hal_size_t map(const DefaultSegmentIterator* source,
                         std::vector<MappedSegmentConstPtr>& results,
                         const Genome* tgtGenome,
                         const std::set<const Genome*>* genomesOnPath,
                         bool doDupes);

protected:
   friend class counted_ptr<DefaultMappedSegment>;
   friend class counted_ptr<const DefaultMappedSegment>;

protected:

   DefaultMappedSegment(SegmentIteratorConstPtr source,
                        SegmentIteratorConstPtr target);

   
   static 
   hal_size_t mapRecursive(const Genome* prevGenome,
                           std::vector<DefaultMappedSegmentConstPtr>& input,
                           std::vector<DefaultMappedSegmentConstPtr>& results,
                           const Genome* tgtGenome,
                           const std::set<const Genome*>* genomesOnPath,
                           bool doDupes);
   static 
   hal_size_t mapUp(DefaultMappedSegmentConstPtr mappedSeg, 
                    std::vector<DefaultMappedSegmentConstPtr>& results);
   static 
   hal_size_t mapDown(DefaultMappedSegmentConstPtr mappedSeg, 
                      hal_size_t childIndex,
                      std::vector<DefaultMappedSegmentConstPtr>& results);
   static 
   hal_size_t mapSelf(DefaultMappedSegmentConstPtr mappedSeg, 
                      std::vector<DefaultMappedSegmentConstPtr>& results);
   
   TopSegmentIteratorConstPtr targetAsTop() const;
   BottomSegmentIteratorConstPtr targetAsBottom() const;
   TopSegmentIteratorConstPtr sourceAsTop() const;
   BottomSegmentIteratorConstPtr sourceAsBottom() const;
   SegmentIteratorConstPtr sourceCopy() const;
   
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

}
#endif
