/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALPOSITIONCACHE_H
#define _HALPOSITIONCACHE_H

#include <map>
#include <string>
#include <vector>
#include <cassert>
#include "halDefs.h"

namespace hal {

/** keep track of bases by storing 2d intervals 
 * For example, if we want to flag positions in a genome
 * that we have visited, this structure will be fairly 
 * efficient provided positions are clustered into intervals */
class PositionCache
{
public:
   PositionCache() : _size(0), _prev(_set.begin()) {}
   PositionCache(const PositionCache &positionCache) : _set(*positionCache.getIntervalSet()), _size(positionCache.size()), _prev(_set.begin()) {}
  // sorted by last index, so each interval is (last, first)
   typedef std::map<hal_index_t, hal_index_t> IntervalSet;
 
   bool insert(hal_index_t pos);
   bool find(hal_index_t pos) const;
   void clear();
   bool check() const;
   hal_size_t size() const { return _size; }
   hal_size_t numIntervals() const { return _set.size(); }

   const IntervalSet* getIntervalSet() const { return &_set; }

protected:

   IntervalSet _set;
   hal_size_t _size;
   IntervalSet::iterator _prev;
};

}

#endif

