/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALAVERAGE_H
#define _HALAVERAGE_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include "hal.h"

namespace hal {

/** Convenience class to keep track of very simple summary stats
 */
template <typename T>
class Average
{
public:
   Average();
   virtual ~Average();
   void reset();
   T getMin() const;
   T getMax() const;
   T getSum() const;
   double getMean() const;
   hal_size_t getCount() const;
   void add(T val, hal_size_t count = 1);   
   Average& operator+=(const Average& other);
   Average& operator/=(hal_size_t N);

protected:
   T _min;
   T _max;
   T _sum;
   hal_size_t _count;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Average<T>& avg);

template <typename T>
inline Average<T>::Average()
{
  reset();
}

template <typename T>
inline Average<T>::~Average()
{
}

template <typename T>
inline T Average<T>::getMin() const
{
  return _min;
}

template <typename T>
inline T Average<T>::getMax() const
{
  return _max;
}

template <typename T>
inline T Average<T>::getSum() const
{
  return _sum;
}

template <typename T>
inline double Average<T>::getMean() const
{
  return (double)_sum / (double)_count;
}

template <typename T>
inline hal_size_t Average<T>::getCount() const
{
  return _count;
}

template <typename T>
inline void Average<T>::add(T val, hal_size_t count)
{
  assert(val > 0);
  _sum += val;
  _count += count;
  _min = std::min(_min, val);
  _max = std::max(_max, val);
}

template <typename T>
inline void Average<T>::reset()
{
  _sum = 0;
  _count = 0;
  _min = std::numeric_limits<T>::max();
  _max = std::numeric_limits<T>::min();
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Average<T>& avg)
{
  os << avg.getMin() << ", " << avg.getMax() << ", " << avg.getMean();
  return os;
}

template <typename T>
inline Average<T>& Average<T>::operator+=(const Average<T>& other)
{
  _min += other._min;
  _max += other._max;
  _sum += other._sum;
  _count += other._count;
  return *this;
}

template <typename T>
inline Average<T>& Average<T>::operator/=(hal_size_t N)
{
  _min /= (T)N;
  _max /= (T)N;
  _sum /= (T)N;
  _count /= (T)N;
  return *this;
}

}

#endif
