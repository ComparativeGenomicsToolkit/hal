/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTYPES_H
#define _HALTYPES_H

namespace hal {

/** 
 *  Keep nontrivial basic types typedefed in this file. 
 */

/**
 * An index in a genome (must be large enough to represent the
 * number of bases in the largest genome in the file)
 */
typedef int64_t genidx_t;

/**
 * An index in a segment (top or bottom) array
 */
typedef int64_t segidx_t;

}
#endif
