/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDNAITERATOR_H
#define _HALDNAITERATOR_H

#include "halDefs.h"
#include "halGenome.h"

namespace hal {

/** 
 * Interface for general dna iterator
 */
class DNAIterator
{
public:
   virtual ~DNAIterator() = 0;
   virtual hal_dna_t getChar() const = 0;
   virtual std::string getString() const = 0;
   virtual void next() const = 0;
   virtual void prev() const = 0;
   virtual void jumpTo(hal_size_t index) const = 0;
   virtual void setChar(hal_dna_t c) = 0;
   virtual void setString(std::string* s) = 0;
   virtual GenomePtr getGenome() = 0;
   virtual GenomeConstPtr getGenome() const = 0;
};

}
#endif
