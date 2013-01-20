/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */


/** Based off of (what I assume to be free code) from 
 * http://ootips.org/yonat/4dev/smart-pointers.html
 * http://ootips.org/yonat/4dev/counted_ptr.h
 *
 * but modified to support static casting and non-const to const casting.  
 * So we should  currently allow assigning 
 * from counted_ptr<X*> to counted_ptr<const X*>, as well as from
 * counted_ptr<DERIVED*> to counted_ptr<BASE*> (and const variations
 * therof).  Dynamic (ie from BASE to DERIVED) casting could be added in
 * later but it would require a new function like counted_ptr.downcast(..).
 *
 * Standardized reference counting smart pointers can be found in Boost and
 * newer gcc's but I do not want to introduce a dependency on either of these
 * for such a small issue.  
 *
 *
 * counted_ptr - simple reference counted pointer.
 *
 * The is a non-intrusive implementation that allocates an additional
 * int and pointer for every counted object.
 */

#ifndef COUNTED_PTR_H
#define COUNTED_PTR_H

namespace hal {


template<typename X>
struct ptr_counter {
   ptr_counter(X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
   X*          ptr;
   unsigned    count;
}; 

template <class X> 
class counted_ptr
{
public:
   template <typename XX> friend class counted_ptr;
   typedef ptr_counter<X> counter;

   explicit counted_ptr(X* p = 0) // allocate a new counter
     : itsCounter(0) {
     if (p) itsCounter = new counter(p);
   }

   ~counted_ptr() {
     release();
   }

   // important to keep the default copy constructor overloaded
   // (even if the templated version below seems to cover it) 
   // to prevent the compiler from making its own and breaking everything
   counted_ptr(const counted_ptr& r) {
     acquire<X>(r.itsCounter);
   }
   
   // important to keep the default assignment operator overloaded
   // (even if the templated version below seems to cover it) 
   // to prevent the compiler from making its own and breaking everything
   counted_ptr& operator=(const counted_ptr& r) {
     if (this != &r)
     {
       release();
       acquire<X>(r.itsCounter);
     }
     return *this;
   }
   
   template <class Y> 
   counted_ptr(const counted_ptr<Y>& r) {
     acquire<Y>(r.itsCounter);
   }

   template <class Y> 
   counted_ptr& operator=(const counted_ptr<Y>& r) {
     if ((const void*)this != (const void*)(&r))
     {
       release();
       acquire<Y>(r.itsCounter);
     }
     return *this;
   }
  
   X& operator*()  const   {return *itsCounter->ptr;}
   X* operator->() const   {return itsCounter->ptr;}
   X* get()        const   {return itsCounter ? itsCounter->ptr : 0;}
   bool unique()   const {
     return (itsCounter ? itsCounter->count == 1 : true);
   }

private:

   counter* itsCounter;

   template <typename Y>
   void acquire(ptr_counter<Y>* c) {
     if (c) {
       (void)static_cast<X*>(c->ptr);
       itsCounter = reinterpret_cast<counter*>(c);
       ++c->count;
     }
     else {
       itsCounter = NULL;
     }
   }

   void release() {
     if (itsCounter) {
       if (--itsCounter->count == 0) {
         delete itsCounter->ptr;
         delete itsCounter;
       }
       itsCounter = 0;
     }
   }
};

// need to specialize on constant template parameters so we can 
// properly account for the const to non-const assignments.  this is
// the only way i could figure out how to do it anyway, since the 
// regular template just seems to understand types and not consts...

// trick to remove const from type (ie template typename).
template <typename T>
struct hal_remove_const {
   typedef T type;
};
template <typename T>
struct hal_remove_const<const T> {
   typedef T type;
};

template <class X> 
class counted_ptr<const X>
{
public:
   template <typename XX> friend class counted_ptr;
   typedef typename hal_remove_const<const X>::type Xnc;
   typedef ptr_counter<Xnc> counter;

   explicit counted_ptr(const X* p = 0) // allocate a new counter
     : itsCounter(0) {
     if (p) itsCounter = new counter(const_cast<Xnc*>(p));
   }

   ~counted_ptr() {
     release();
   }

   // important to keep the default copy constructor overloaded
   // (even if the templated version below seems to cover it) 
   // to prevent the compiler from making its own and breaking everything
   counted_ptr(const counted_ptr<const X>& r) {
     acquire<Xnc>(r.itsCounter);
   }
   counted_ptr(const counted_ptr<X>& r) {
     acquire<Xnc>(const_cast<counter*>(r.itsCounter));
   }
   
   // important to keep the default assignment operator overloaded
   // (even if the templated version below seems to cover it) 
   // to prevent the compiler from making its own and breaking everything
   counted_ptr<const X>& operator=(const counted_ptr<const X>& r) {
     if (this != &r)
     {
       release();
       acquire<Xnc>(r.itsCounter);
     }
     return *this;
   }
   counted_ptr<const X>& operator=(
     const counted_ptr<X>& r) {
     return this->operator=(const_cast<const counted_ptr<const X>&>(r));
   }
  
   template <class Y> 
   counted_ptr(const counted_ptr<Y>& r) {
     acquire<Y>(r.itsCounter);
   }

   template <class Y> 
   counted_ptr& operator=(const counted_ptr<Y>& r) {
     if ((const void*)this != (const void*)(&r))
     {
       release();
       acquire<Y>(r.itsCounter);
     }
     return *this;
   }
   
   const X& operator*()  const   {return *itsCounter->ptr;}
   const X* operator->() const   {return itsCounter->ptr;}
   const X* get()        const   {return itsCounter ? itsCounter->ptr : 0;}
   bool unique()   const {
     return (itsCounter ? itsCounter->count == 1 : true);
   }

private:

   counter* itsCounter;

   template <typename Y>
   void acquire(ptr_counter<Y>* c) {
     if (c) {
       (void)static_cast<const X*>(c->ptr);
       itsCounter = reinterpret_cast<counter*>(c);
       ++c->count;
     }
     else {
       itsCounter = NULL;
     }
   }

   void release() {
     if (itsCounter) {
       if (--itsCounter->count == 0) {
         delete itsCounter->ptr;
         delete itsCounter;
       }
       itsCounter = 0;
     }
   }
};


}
#endif // COUNTED_PTR_H
