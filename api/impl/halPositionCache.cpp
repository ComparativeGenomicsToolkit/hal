/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <sstream>
#include "halPositionCache.h"
#include "hal.h"

using namespace std;
using namespace hal;


bool PositionCache::insert(hal_index_t pos)
{
  IntervalSet::iterator i = _set.lower_bound(pos);
  IntervalSet::iterator j;
  if (i != _set.end() && i->second <= pos)
  {
    assert(i->first >= pos);
    return false;
  }

  //merge to beginning of existing interval
  if (i != _set.end() && i->second == pos + 1)
  {
    --i->second;
  }
  else
  {
    // create new unit interval
    i = _set.insert(pair<hal_index_t, hal_index_t>(pos, pos)).first;
  }
  assert(i->second <= i->first);
  // merge abutting left interval
  if (i != _set.begin())
  {
    j = i;
    --j;
    if (j->first == i->second - 1)
    {
      i->second = j->second;
      assert(i->second <= i->first);
      _set.erase(j);
    }    
  }

  // merge abutting right interval
  j = i;
  ++j;
  if ( j != _set.end() && j->second == i->first + 1)
  {
     j->second = i->second;
     assert(j->second <= j->first);
     _set.erase(i);
  }

  ++_size;
  assert(find(pos) == true);
  return true;
}

bool PositionCache::find(hal_index_t pos) const
{
  IntervalSet::const_iterator i = _set.lower_bound(pos);
  if (i != _set.end() && i->second <= pos)
  {
    return true;
  }
  return false;
}

void PositionCache::clear()
{
  _set.clear();
  _size = 0;
}

// for debugging
bool PositionCache::check() const
{
  hal_size_t size = 0;
  for (IntervalSet::const_iterator i = _set.begin(); i != _set.end(); ++i)
  {
    size += (i->first + 1) - i->second;
    IntervalSet::const_iterator j = i;
    ++j;
    if (j != _set.end())
    { 
      // test overlap
      if (j->second <= i->first || i->second >= j->first)
      {
        return false;
      }
      // test merge
      if (j->second == i->first + 1)
      {
        return false;
      }
    }
  }
  return size == _size;
}
