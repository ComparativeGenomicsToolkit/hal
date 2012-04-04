/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGENOME_H
#define _HALGENOME_H

#include "halDefs.h"
#include "halAlignment.h"
#include "halSegmentIterator.h"
#include "halDNAIterator.h"

namespace hal {

/** 
 * Abstract base class for a genome
 */
class Genome
{
public:
   virtual ~Genome() = 0;
     
   virtual const std::string& getName() const = 0;
   virtual AlignmentPtr getAlignment() = 0;
   virtual AlignmentConstPtr getAlignment() const = 0;
   virtual hal_size_t getSequenceLength() const = 0;
   virtual hal_size_t getNumberTopSegments() const = 0;
   virtual hal_size_t getNumberBottomSegmentes() const = 0;
   virtual SegmentIteratorPtr getSegmentIterator(hal_bool_t top, 
                                              hal_index_t position) = 0;
   virtual SegmentIteratorConstPtr getSegmentIterator(
     hal_bool_t top, hal_index_t position) const = 0;
   
   virtual MetaDataPtr getMetaData() = 0;
   virtual MetaDataConstPtr getMetaData() const = 0;
  
};

}
#endif
