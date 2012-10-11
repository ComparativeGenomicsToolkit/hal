/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEQUENCEITERATOR_H
#define _HALSEQUENCEITERATOR_H

#include "halDefs.h"
#include "halSequence.h"

namespace hal {

/** 
 * Iterate over sequences in the genome.  
 */
class SequenceIterator
{
public:
   /** Create a duplicate iterator referring to the same sequence
    * which itself is not copied */
   virtual SequenceIteratorPtr copy() = 0;

   /** Create a duplicate iterator referring to the same sequence
    * which itself is not copied */
   virtual SequenceIteratorConstPtr copy() const = 0;

   /** Move iterator to next sequence in the genome */
   virtual void toNext() const = 0;

   /** Move iterator to previous sequence in the genome */
   virtual void toPrev() const = 0;
   
   /** Return pointer to the sequence */
   virtual Sequence* getSequence() = 0;

   /** Return pointer to the sequence */
   virtual const Sequence* getSequence() const = 0;

   /** Test if iterator points to same sequence as other iterator */
   virtual bool equals(SequenceIteratorConstPtr p2) const = 0;

protected:
   friend class counted_ptr<SequenceIterator>;
   friend class counted_ptr<const SequenceIterator>;
   virtual ~SequenceIterator() = 0;
};

inline SequenceIterator::~SequenceIterator() {}

inline bool operator==(SequenceIteratorConstPtr p1,
                       SequenceIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(SequenceIteratorConstPtr p1,
                       SequenceIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}
#endif
