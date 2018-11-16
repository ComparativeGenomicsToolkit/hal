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
   /** Destructor */
    virtual ~SequenceIterator() {
    }

   /** Create a duplicate iterator referring to the same sequence
    * which itself is not copied */
   virtual SequenceIteratorPtr clone() const = 0;

   /** Move iterator to next sequence in the genome */
   virtual void toNext() = 0;

   /** Move iterator to previous sequence in the genome */
   virtual void toPrev() = 0;

    /** has the iterator reach the end of the traversal in the direction of
     * movement? */
    virtual bool atEnd() const = 0;
    
   /** Return pointer to the sequence */
   virtual const Sequence* getSequence() const = 0;

   /** Test if iterator points to same sequence as other iterator */
   virtual bool equals(SequenceIteratorPtr p2) const = 0;
};

inline bool operator==(SequenceIteratorPtr p1,
                       SequenceIteratorPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(SequenceIteratorPtr p1,
                       SequenceIteratorPtr p2)
{
  return !(p1 == p2);
}

}
#endif
