/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DYNBOTTOMSEGMENT_H
#define _DYNBOTTOMSEGMENT_H

#include "bottomSegment.h"

namespace hal {

/** 
 *  General implementation of a "Bottom" DNA segment
 *  Useful for larger numbers of children when we can't
 *  compile a zillion template instances of staticBottomSegment..
 */
class DynBottomSegment: public BottomSegment
{
public:
	/** Number of descent edges to children (any of these
	 * can still be null values */
	int getNumChildren() const;
	/** Get segment index for top array of ith child array */
	segidx_t getChildIndex(int i) const;
	/** Get reverse complement flag for ith child index */
	bool getReverseComp(int i) const;

protected:
	segidx_t* _childSegIndex;
	bool* _reverseComplement;
};  

}
#endif
