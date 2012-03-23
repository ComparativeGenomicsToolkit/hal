/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _STATICBOTTOMSEGMENT_H
#define _STATICBOTTOMSEGMENT_H

#include "bottomSegment.h"

namespace hal {

/** 
 *  Templated implementation of a "Bottom" DNA segment
 *  Useful for small number of children (saves calls to malloc)
 */
template <int NumChildren>
class StaticBottomSegment: public BottomSegment
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
	segidx_t _childSegIndex[NumChildren];
	bool _reverseComplement[NumChildren];
};  

}
#endif
