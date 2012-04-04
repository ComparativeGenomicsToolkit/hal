/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */


/** Based off of (what I assume to be free code) from 
 * http://ootips.org/yonat/4dev/smart-pointers.html
 * http://ootips.org/yonat/4dev/counted_ptr.h
 *
 * but (poorly) hacked to support const casting.  Should replace with
 * a proper implementation, but using the standard shared_ptr<>
 * is pretty heavy on the dependency front (only in super new c++ 
 * compilers and boost). 
 */

/*
 * counted_ptr - simple reference counted pointer.
 *
 * The is a non-intrusive implementation that allocates an additional
 * int and pointer for every counted object.
 */

#ifndef COUNTED_PTR_H
#define COUNTED_PTR_H

/* For ANSI-challenged compilers, you may want to #define
 * NO_MEMBER_TEMPLATES or explicit */

namespace hal {

template <typename T>
struct counted_ptr_remove_const
{
    typedef T type;
};
template <typename T>
struct counted_ptr_remove_const<const T>
{
    typedef T type;
};

template<typename X>
struct ptr_counter {
   ptr_counter(X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
   X*          ptr;
   unsigned    count;
}; 

template<typename X>
struct ptr_counter<const X> {
   ptr_counter(const X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
   const X*          ptr;
   unsigned    count;
   
   ptr_counter& operator=(const ptr_counter<typename counted_ptr_remove_const<X>::type>& other) {
     ptr = other.ptr;
     count = other.count;
   }
};
template <class X> class counted_ptr
{
public:
   friend class counted_ptr<const X>;
   typedef X element_type;
   typedef ptr_counter<X> counter;

    explicit counted_ptr(X* p = 0) // allocate a new counter
      : itsCounter(0) {if (p) itsCounter = new counter(p);}
    ~counted_ptr()
        {release();}
    counted_ptr(const counted_ptr& r) throw()
        {acquire(r.itsCounter);}
    counted_ptr& operator=(const counted_ptr& r)
    {
        if (this != &r) {
            release();
            acquire(r.itsCounter);
        }
        return *this;
    } 


    /*
#ifndef NO_MEMBER_TEMPLATES
    template <class Y> friend class counted_ptr<Y>;
    template <class Y> counted_ptr(const counted_ptr<Y>& r) throw()
        {acquire(r.itsCounter);}
    template <class Y> counted_ptr& operator=(const counted_ptr<Y>& r)
    {
        if (this != &r) {
            release();
            acquire(r.itsCounter);
        }
        return *this;
    }
#endif // NO_MEMBER_TEMPLATES
    */

    X& operator*()  const throw()   {return *itsCounter->ptr;}
    X* operator->() const throw()   {return itsCounter->ptr;}
    X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
    bool unique()   const throw()
        {return (itsCounter ? itsCounter->count == 1 : true);}

private:

   counter* itsCounter;

    void acquire(counter* c) throw()
    { // increment the count
        itsCounter = c;
        if (c) ++c->count;
    }

    void release()
    { // decrement the count, delete if it is 0
        if (itsCounter) {
            if (--itsCounter->count == 0) {
                delete itsCounter->ptr;
                delete itsCounter;
            }
            itsCounter = 0;
        }
    }
};


template <class X> class counted_ptr<const X>
{
public:
    typedef const X element_type;
   typedef ptr_counter<const X> counter;

    explicit counted_ptr(const X* p = 0) // allocate a new counter
        : itsCounter(0) {if (p) itsCounter = new counter(p);}
   ~counted_ptr()
        {release();}
    counted_ptr(const counted_ptr& r) throw()
        {acquire(r.itsCounter);}
    counted_ptr& operator=(const counted_ptr& r)
    {
        if (this != &r) {
            release();
            acquire(r.itsCounter);
        }
        return *this;
    }

   counted_ptr(const counted_ptr<typename counted_ptr_remove_const<X>::type>& r) throw()
      {acquire(r.itsCounter);}
   counted_ptr& operator=(const counted_ptr<typename counted_ptr_remove_const<X>::type>& r)
   {
     if (this != &r) {
       release();
       acquire(r.itsCounter);
     }
     return *this;
   }
/*
#ifndef NO_MEMBER_TEMPLATES
    template <class Y> friend class counted_ptr<Y>;
    template <class Y> counted_ptr(const counted_ptr<Y>& r) throw()
        {acquire(r.itsCounter);}
    template <class Y> counted_ptr& operator=(const counted_ptr<Y>& r)
    {
        if (this != &r) {
            release();
            acquire(r.itsCounter);
        }
        return *this;
    }
#endif // NO_MEMBER_TEMPLATES
*/
    const X& operator*()  const throw()   {return *itsCounter->ptr;}
    const X* operator->() const throw()   {return itsCounter->ptr;}
    const X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
    bool unique()   const throw()
        {return (itsCounter ? itsCounter->count == 1 : true);}

private:

   counter* itsCounter;

    void acquire(counter* c) throw()
    { // increment the count
        itsCounter = c;
        if (c) ++c->count;
    }

    void acquire(ptr_counter<typename counted_ptr_remove_const<X>::type>* c) throw()
    { // increment the count
      itsCounter = (ptr_counter<const X>*)c;
        if (c) ++c->count;
    }

    void release()
    { // decrement the count, delete if it is 0
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
