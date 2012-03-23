/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _BOTTOMSEGMENT_H
#define _BOTTOMSEGMENT_H

#include "genomeSegment.h"

namespace hal {

/** 
 *  Abstract base class for a "Bottom" DNA segment
 */
class BottomSegment: public GenomeSegment<genidx_t, segidx_t>
{
public:
	/** Number of descent edges to children (any of these
	 * can still be null values */
	virtual int getNumChildren() const = 0;
	/** Get segment index for top array of ith child array */
	virtual segidx_t getChildIndex(int i) const = 0;
	/** Get reverse complement flag for ith child index */
	virtual bool getReverseComp(int i) const = 0;
	/** Get next repeat segment's index */
	segidx_t getRepeatIndex() const;

protected:
	segidx_t _repeatIndex;
};  

}
#endif
