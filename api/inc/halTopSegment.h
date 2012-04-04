/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENT_H
#define _HALTOPSEGMENT_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for a top segment of DNA
 */
class TopSegment : public Segment
{
public:
  virtual ~TopSegment() = 0;
  virtual hal_index_t getParentIndex() const = 0;
  virtual void setParentIndex(hal_index_t parIdx) = 0;
  virtual hal_bool_t getParentReversed() const = 0;
  virtual void setParentReversed(hal_bool_t isReversed) = 0;
  virtual hal_index_t getBottomParseIndex() const = 0;
  virtual void setBottomParseIndex(hal_index_t botParseIdx) = 0;
  virtual hal_offset_t getBottomParseOffset() const = 0;
  virtual void setBottomParseOffset(hal_offset_t botParseOffset) = 0;
};


}
#endif
