/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDNAITERATOR_H
#define _HALDNAITERATOR_H

#include "halDefs.h"

namespace hal {

/** 
 * Interface for general dna iterator
 */
class DNAIterator
{
public:

   /** Get the DNA character at this position (if revsersed is set
    * the reverse compelement is returned */
   virtual char getChar() const = 0;

   /** Set the DNA character at this position (if revsersed is set
    * the reverse compelement is stored
    * @param c DNA character to set */
   virtual void setChar(char c) = 0;

   /** Move to previous position (equiv. to toRight if reversed)*/
   virtual void toLeft() const = 0;

   /** Move to next position (equiv. to toLeft if reversed)*/
   virtual void toRight() const = 0;

   /** Jump to any point on the genome (can lead to 
    * inefficient paging from disk if used irresponsibly)
    * @param index position in array to jump to */
   virtual void jumpTo(hal_size_t index) const = 0;

   /** switch to base's reverse complement */
   virtual void toReverse() const = 0;

   /** Check whether iterator is on base's complement */
   virtual bool getReversed() const = 0;

   /** Set the iterator's reverse complement status */
   virtual void setReversed(bool reversed) const = 0;
     
   /** Get the containing (read-only) genome */
   virtual const Genome* getGenome() const = 0;

   /** Get the containing genome */
   virtual Genome* getGenome() = 0;

   /** Get the containing (read-only) sequence */
   virtual const Sequence* getSequence() const = 0;

   /** Get the containing sequence */
   virtual Sequence* getSequence() = 0;

   /** Get the index of the base in the dna array */
   virtual hal_index_t getArrayIndex() const = 0;

   /** Compare (array indexes) of two iterators */
   virtual bool equals(DNAIteratorConstPtr& other) const = 0;

   /** Compare (array indexes) of two iterators */
   virtual bool leftOf(DNAIteratorConstPtr& other) const = 0;

protected:

   friend class counted_ptr<DNAIterator>;
   friend class counted_ptr<const DNAIterator>;
   virtual ~DNAIterator() = 0;
};

inline DNAIterator::~DNAIterator() {}

}

#ifndef NDEBUG
#include "halGenome.h"
namespace hal {
inline std::ostream& operator<<(std::ostream& os, 
                                const DNAIterator* dna)
{
  const Genome* genome = dna->getGenome();
  os << "dna: ";
  os << "Genome=" << genome->getName();
  os << " Seq=" << dna->getSequence()->getName();
  os << " idx=" << dna->getArrayIndex();

  if (dna->getArrayIndex() >= 0 && 
      dna->getArrayIndex() < (hal_index_t)genome->getSequenceLength())
  {
    os << " val=" << dna->getChar();
  }
  return os;
}
}
#endif

#endif
