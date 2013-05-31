/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAPPEDSEGMENT_H
#define _HALMAPPEDSEGMENT_H

#include "halDefs.h"
#include "halSlicedSegment.h"

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

   /** Get the original segment from which this segment was mapped */
   virtual SlicedSegmentConstPtr getSource() const = 0;

   /** Comparison used to store in stl sets and maps.  We sort based
    * on the coordinate of the mapped segemnt's (target) interval as the primary
    * index and the target genome as the secondary index.  */
   virtual bool lessThan(const MappedSegmentConstPtr& other) const = 0;

   /** Comparison used to store in STL containers.  We sort based
    * on the coordinate of the *Source* genome interval as the primary
    * index and the target genome as the secondary index.  */
   virtual bool lessThanBySource(const MappedSegmentConstPtr& other) const = 0;

   /** Comparison used to determine uniqueness in lists.  Tests lessThan
    * in both directions.  Note that equality is the same regardless of 
    * whether or not we use the source segment as our primary index. */
   virtual bool equals(const MappedSegmentConstPtr& other) const = 0;

   /** Flip the mapping direction.  This segment becomes the source, and
    * the source becomes this.*/
   virtual void flip() const = 0;

   /** Reverse both segments.  Also swap their start and end offsets. 
    * Note that toReverse() does not reverse the source segment.*/
   virtual void fullReverse() const = 0;
   
   /** Return of a copy of the mapped segment */
   virtual MappedSegmentConstPtr copy() const = 0;

   /** Test if mapped segment can be merged to the right with input 
    * segment.  will return false if the right coordinate of this is in
    * either (optional) cutSet.*/
   virtual bool canMergeRightWith(
     const MappedSegmentConstPtr& next,
     const std::set<hal_index_t>* cutSet = NULL,
     const std::set<hal_index_t>* sourceCutSet = NULL) const = 0;

   /** Functor for sorted STL containers, sorting by origin as primary 
    * index */
   struct LessSource { bool operator()(const MappedSegmentConstPtr& ms1,
                                       const MappedSegmentConstPtr& ms2) const {
     return ms1->lessThanBySource(ms2); }
   };

   /** Functor for sorted STL containers, sorting  by target as primary 
    * index */
   struct Less { bool operator()(const MappedSegmentConstPtr& ms1,
                                 const MappedSegmentConstPtr& ms2) const {
     return ms1->lessThan(ms2); }
   };
 
   /** Functor for STL sorted lists to test for uniqueness */
   struct EqualTo { bool operator()(const MappedSegmentConstPtr& ms1,
                                    const MappedSegmentConstPtr& ms2) const {
     return ms1->equals(ms2); }
   };

   // NEEDS TO BE ADDED TO SEGMENT INTERFACE
   virtual void print(std::ostream& os) const = 0;

protected:
   friend class counted_ptr<MappedSegment>;
   friend class counted_ptr<const MappedSegment>;
   virtual ~MappedSegment() = 0;
};

inline MappedSegment::~MappedSegment() {}

}

namespace std {
template<>
struct less<hal::MappedSegmentConstPtr>
{
   bool operator()(const hal::MappedSegmentConstPtr& m1,
                   const hal::MappedSegmentConstPtr& m2) const
      {
        return m1->lessThan(m2);
      }
};

}

#endif
