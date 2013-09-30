/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef COUNTED_PTR_H
#define COUNTED_PTR_H

namespace hal {

// trick to remove const from type (ie template typename).
// (we keep non-const pointers internally to simplify some thigns)
template <typename T>
struct hal_remove_const {
   typedef T type;
};
template <typename T>
struct hal_remove_const<const T> {
   typedef T type;
};

/** 
 * Standardized reference counting smart pointers can be found in Boost and
 * newer gcc's but I do not want to introduce a dependency on either of these
 * for such a small issue.  hence the hack below:
 *
 * Casting should work sensibly.  Const -> non-const and derived -> base
 * supported automatically.  Dynamic casting (base to derived) 
 * is achieved with the downCast method().
 *
 */
template <class T> 
class counted_ptr
{
public:
   template <typename U> friend class counted_ptr;
   typedef typename hal_remove_const<const T>::type Tnc;

   /** Allocate a new counter first reference to pointer */
   explicit counted_ptr(T* p = 0); 

   ~counted_ptr();

   /**important to keep the default copy constructor overloaded
    * (even if the templated version below seems to cover it) 
    * to prevent the compiler from making its own and breaking everything */
   counted_ptr(const counted_ptr<T>& r);
   
   /**important to keep the default assignment operator overloaded
    * (even if the templated version below seems to cover it) 
    * to prevent the compiler from making its own and breaking everything */
   counted_ptr& operator=(const counted_ptr<T>& r);
   
   template <class U> 
   counted_ptr(const counted_ptr<U>& r);

   template <class U> 
   counted_ptr& operator=(const counted_ptr<U>& r);
  
   T& operator*() const;
   T* operator->() const;
   T* get() const;
   bool unique() const;

   /** Should be equivalent to a dynami_cast */
   template <typename U>
   U downCast();

private:

   /** Absorb the pointer and counter from another smart pointer. 
    * Counter gets increased */
   template <typename U>
   void acquire(const counted_ptr<U>& c);

   /** As above but tries to dynamic_cast.  Unlike above, there is no 
    * compile-time check as to whether it will works.  If the cast is
    * incompatible, then NULL is returned */
   template <typename U>
   void dynamicAcquire(const counted_ptr<U>& c);

   void release();

private:

   unsigned* _counter;
   Tnc* _ptr;
};



template <class T> 
inline counted_ptr<T>::counted_ptr(T* p) : _ptr(const_cast<Tnc*>(p)) 
{
  _counter = p ? new unsigned(1) : NULL;
}

template <class T> 
inline counted_ptr<T>::~counted_ptr() 
{
  release();
}

template <class T> 
inline counted_ptr<T>::counted_ptr(const counted_ptr<T>& r) 
{
  acquire<T>(r);
}
   
template <class T> 
inline counted_ptr<T>& counted_ptr<T>::operator=(const counted_ptr<T>& r) 
{
  if (this != &r)
  {
    release();
    acquire<T>(r);
  }
  return *this;
}
   
template <class T> 
template <class U>
inline counted_ptr<T>::counted_ptr(const counted_ptr<U>& r) 
{
  acquire<U>(r);
}

template <class T> 
template <class U>
inline counted_ptr<T>& counted_ptr<T>::operator=(const counted_ptr<U>& r) 
{
  if ((const void*)this != (const void*)(&r))
  {
    release();
    acquire<U>(r);
  }
  return *this;
}
  
template <class T> 
inline T& counted_ptr<T>::operator*() const
{
  return *_ptr;
}

template <class T> 
inline T* counted_ptr<T>::operator->() const
{
  return _ptr;
}

template <class T> 
inline T* counted_ptr<T>::get() const
{
  return _ptr;
}

template <class T> 
inline bool counted_ptr<T>::unique() const 
{
  return (_ptr ? *_counter == 1 : true);
}

template <class T> 
template <typename U>
inline U counted_ptr<T>::downCast() 
{
  U other;
  other.dynamicAcquire(*this);
  return other;
}

template <class T> 
template <typename U>
inline void counted_ptr<T>::acquire(const counted_ptr<U>& c) 
{
  if (c._ptr) 
  {
    _ptr = const_cast<Tnc*>(static_cast<T*>(c._ptr));
    _counter = c._counter;
    ++(*_counter);
  }
  else 
  {
    _counter = NULL;
    _ptr = NULL;
  }
}


template <class T> 
template <typename U>
inline void counted_ptr<T>::dynamicAcquire(const counted_ptr<U>& c) 
{
  Tnc* temp = const_cast<Tnc*>(dynamic_cast<T*>(c._ptr));
  if (temp) 
  {
    _ptr = temp;
    _counter = c._counter;
    ++(*_counter);
  }
  else 
  {
    _counter = NULL;
    _ptr = NULL;
  }
}

template <class T> 
inline void counted_ptr<T>::release() 
{
  if (_counter) 
  {
    if (--(*_counter) == 0) 
    {
      delete _ptr;
      delete _counter;
    }
    _ptr = 0;
    _counter = 0;
  }
}

template <class T>
inline std::ostream& operator<<(std::ostream& os, const counted_ptr<T>& cp)
{
  return os << cp.get();
}

}
#endif // COUNTEDPTR_H
