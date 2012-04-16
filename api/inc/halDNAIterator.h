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

   /** Get the DNA character at this position */
   virtual hal_dna_t getChar() const = 0;

   /** Get the reverse Complemenet of the DNA character at this position */
   virtual hal_dna_t getCompChar() const = 0;

   /** Set the DNA character at this position 
    * @param c DNA character to set */
   virtual void setChar(hal_dna_t c) = 0;

   /** Move to previous position */
   virtual void toLeft() const = 0;

   /** Move to next position */
   virtual void toRight() const = 0;

   /** Jump to any point on the genome (can lead to 
    * inefficient paging from disk if used irresponsibly)
    * @param index position in array to jump to */
   virtual void jumpTo(hal_size_t index) const = 0;
     
   /** Get the containing (read-only) genome */
   virtual const Genome* getGenome() const = 0;

   /** Get the containing genome */
   virtual Genome* getGenome() = 0;

   /** Get the index of the segment in the segment array */
   virtual hal_index_t getArrayIndex() const = 0;

protected:

   friend class counted_ptr<DNAIterator>;
   friend class counted_ptr<const DNAIterator>;
   virtual ~DNAIterator() = 0;
};

inline DNAIterator::~DNAIterator() {}
}
#endif
