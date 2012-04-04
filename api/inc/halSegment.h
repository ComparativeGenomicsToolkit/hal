/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENT_H
#define _HALSEGMENT_H

#include <string>
#include "halDefs.h"
#include "halGenome.h"

namespace hal {

/** 
 * Interface for a segment of DNA
 */
class Segment
{
public:
   virtual ~Segment() = 0;
   virtual hal_size_t getLength() const = 0;
   virtual void setLength() = 0;
   virtual GenomeConstPtr getGenome() const = 0;
   virtual GenomePtr getGenome() = 0;
   virtual hal_index_t getStartPostion() const = 0;
   virtual void setStartPosition() = 0;
   virtual std::string getSequence() const = 0;
   virtual hal_index_t getNextParalogyIndex() const = 0;
};


}
#endif
