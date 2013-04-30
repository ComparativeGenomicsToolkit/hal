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


class DefaultMappedSegment : virtual public DefaultSegmentIterator,
                             virtual public MappedSegment
{
public:
   
   virtual ~DefaultMappedSegment();
   
   // MAPPED SEGMENT INTERFACE 
   virtual SlicedSegmentConstPtr getSource() const;

   // INTERNAL METHODS
   static hal_size_t map(const DefaultSegmentIterator* source,
                         std::vector<MappedSegmentConstPtr>& results,
                         const Genome* tgtGenome,
                         const std::set<const Genome*>* genomesOnPath,
                         bool doDupes);

protected:

   DefaultMappedSegment(DefaultSegmentIteratorConstPtr source,
                        DefaultSegmentIteratorConstPtr target);

   bool mapUp(std::vector<MappedSegmentConstPtr>& results);
   bool mapDown(std::vector<MappedSegmentConstPtr>& results);

   virtual SegmentPtr getSegment();
   virtual SegmentConstPtr getSegment() const;
   virtual hal_size_t getNumSegmentsInGenome() const;
   
protected:

   SegmentIteratorConstPtr _source;
   SegmentConstPtr _segment;
};

}
#endif
