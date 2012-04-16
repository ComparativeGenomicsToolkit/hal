/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDEFS_H
#define _HALDEFS_H

#include <string>
#include <vector>
#include <stdexcept>
#include "counted_ptr.h"

namespace hal {

/** 
 *  Keep simple compiler-related definitions in this file. 
 */

/**
 * An index in a genome (must be large enough to represent the
 * number of bases in the largest genome in the file)
 */
typedef int64_t hal_index_t;
extern const hal_index_t NULL_INDEX;

/**
 * An index in a segment (top or bottom) array
 */
typedef int64_t hal_offset_t;

/**
 * Size type
 */
typedef uint64_t hal_size_t;

/**
 * A DNA character
 */
typedef char hal_dna_t;

/**
 * Boolean
 */
typedef bool hal_bool_t;

/*
 * General usage exception class, used for all critical errors. 
 */
typedef std::runtime_error hal_exception;

/**
 * Smart pointer type.  Should be official one like boost::shared_ptr
 * but we use this guy's to avoid the dependency for now.  Goes in
 * a struct because typedef doesn't support template arguments otherwise. 
 */
template <typename T>
struct smart_ptr {
   typedef counted_ptr<T> type;
};

// FORWARD DECLARATIONS
#define HAL_FORWARD_DEC_CLASS(T) \
  class T;\
  typedef smart_ptr<T>::type T ## Ptr;\
  typedef smart_ptr<const T>::type T ## ConstPtr;

HAL_FORWARD_DEC_CLASS(Alignment)
HAL_FORWARD_DEC_CLASS(Genome)
HAL_FORWARD_DEC_CLASS(MetaData)
HAL_FORWARD_DEC_CLASS(TopSegment)
HAL_FORWARD_DEC_CLASS(BottomSegment)
HAL_FORWARD_DEC_CLASS(Segment)
HAL_FORWARD_DEC_CLASS(Sequence)
HAL_FORWARD_DEC_CLASS(SegmentIterator)
HAL_FORWARD_DEC_CLASS(TopSegmentIterator)
HAL_FORWARD_DEC_CLASS(BottomSegmentIterator)
HAL_FORWARD_DEC_CLASS(DNAIterator)
HAL_FORWARD_DEC_CLASS(SequenceIterator)
  
}
#endif
