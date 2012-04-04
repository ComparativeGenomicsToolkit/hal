/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENT_H
#define _HALBOTTOMSEGMENT_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for a bottom segment of DNA
 */
class BottomSegment : public Segment
{
public:
  virtual ~BottomSegment() = 0;
  virtual hal_size_t getNumChildren() const = 0;
  virtual hal_index_t getChildIndex(hal_size_t i) const = 0;
  virtual void setChildIndex(hal_size_t child, hal_index_t childIndex) = 0;
  virtual hal_bool_t getChildReversed(hal_size_t i) const = 0;
  virtual void setChildReversed(hal_size_t child, hal_bool_t isReversed) = 0;
  virtual hal_index_t getParentParseIndex() const = 0;
  virtual void setParentParseIndex(hal_index_t parParseIndex) = 0;
  virtual hal_offset_t getParentParseOffset() const = 0;
  virtual void setParentParseOffset(hal_offset_t parParseOffset) = 0;
};


}
#endif
