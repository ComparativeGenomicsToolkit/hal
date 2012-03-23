/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _GENOMESEGMENT_H
#define _GENOMESEGMENT_H

#include "haltypes.h"

namespace hal {

/** 
 * Abstract base class for a DNA segment (interval in genome)
 * Segments for now are glorified structs that hold the coordinates
 * and little else.  More useful interface will be at higher level
 */
class GenomeSegment
{
public:
/** Get start coordinate in genome */
	genidx_t getStart() const;
/** Get length of segment */
	genidx_t getLength() const;
/** Get index of leftmost segment overlapping in same genome */
	segidx_t getParseIdx() const;
/** Get offset of above */
	genidx_t getParseOffset() const;

/** Null index value */
	static const segidx_t NullIndex;

protected:
	genidx_t _start;
	genidx_t _length;
	segidx_t _parseIdx;
	genidx_t _parseOffset;
};  

}
#endif

