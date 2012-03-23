/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _TOPSEGMENT_H
#define _TOPSEGMENT_H

#include "genomeSegment.h"

namespace hal {

/** 
 * "Top" DNA segment.  It will be connected with an ancestral genome
 */
class TopSegment: public GenomeSegment
{
public:
	/** Get segment index in parent genome bottom array */
	segidx_t getParentIndex() const;
	/** Get reverse complement flag for parent index */
	bool getReverseComp() const;
	/** Get next repeat segment's indexf */
	segidx_t getRepeatIndex() const;

protected:
	segidx_t _parentSegIdx;
	bool _parentSegReverse;
	segidx_t _repeatSegIdx;
};  

}
#endif
