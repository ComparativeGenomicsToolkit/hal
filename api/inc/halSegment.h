/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENT_H
#define _HALSEGMENT_H

#include <string>
#include "halDefs.h"

namespace hal {

/** 
 * Interface for a segment of DNA. Note that segments should
 * not be written to outside of creating new genomes.
 */
class Segment
{
public:
 
   /** Get the length of the segment (number of bases) */
   virtual hal_size_t getLength() const = 0;

   /** Set the length of the segment
    * @param length New length of segment */
   virtual void setLength(hal_size_t length) = 0;

   /** Get the containing (read-only) genome */
   virtual const Genome* getGenome() const = 0;

   /** Get the containing genome */
   virtual Genome* getGenome() = 0;

   /** Get the containing (read-only) sequence */
   virtual const Sequence* getSequence() const = 0;

   /** Get the containing sequence */
   virtual Sequence* getSequence() = 0;

   /** Get the segment's start position in the genome */
   virtual hal_index_t getStartPosition() const = 0;

   /** Set the segment's start position in the genome 
    * @param startPos Start position */
   virtual void setStartPosition(hal_index_t startPos) = 0;

   /** Get the index of the segment in the segment array */
   virtual hal_index_t getArrayIndex() const = 0;

   /** Check whether segment is the first segment of a sequence */
   virtual bool isFirst() const = 0;

   /** Check whether segment is the last segment of a sequence */
   virtual bool isLast() const = 0;


protected:
   
   /** Destructor */
   virtual ~Segment() = 0;
};

inline Segment::~Segment() {}

}
#endif
