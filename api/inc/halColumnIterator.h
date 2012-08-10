/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNITERATOR_H
#define _HALCOLUMNITERATOR_H

#include <list>
#include <map>
#include <set>
#include "halDefs.h"
#include "halDNAIterator.h"
#include "halSequence.h"

namespace hal {

/** 
 * Interface Column iterator for allowing traditional maf-like (left-to-right)  
 * parsing of a hal alignment.  Columns are iterated with respect to
 * a spefified reference genome.  This isn't the most efficient way
 * to explore the hal structure, which is designed for bottom-up and/or
 * top-down traversal.  
 */
class ColumnIterator 
{
public:

   // we can compare genomes by pointers (because they are persistent
   // and unique, though it's still hacky) but we can't do the same 
   // for anything else, including sequences.  
   struct SequenceLess { bool operator()(const hal::Sequence* s1,
                                         const hal::Sequence* s2) const {
     return s1->getGenome() < s2->getGenome() || (
       s1->getGenome() == s2->getGenome() && 
       s1->getStartPosition() < s2->getStartPosition()); }
   };

   typedef std::vector<hal::DNAIteratorConstPtr> DNASet;
   typedef std::map<const hal::Sequence*, DNASet*, SequenceLess> ColumnMap;

   /** Move column iterator one column to the right along reference
    * genoem sequence */
   virtual void toRight() const = 0;

   /** Test if column iterator is left  (only in terms
    * of the reference genome) as another iterator.  Used to bound
    * a loop, for example.  in the case of paralogy, the rightmost
    * position of this is compared to the leftmost position of other*/
   virtual bool leftOf(ColumnIteratorConstPtr other) const = 0;
   
   /** Get a pointer to the reference genome for the column iterator */
   virtual const hal::Genome* getReferenceGenome() const = 0;

   /** Get a pointer to the reference sequence for the column iterator */
   virtual const hal::Sequence* getReferenceSequence() const = 0;

   /** Get a pointer to the column map */
   virtual const ColumnMap* getColumnMap() const = 0;

   /** Get the index of the column in the reference genome's array */
   virtual hal_index_t getArrayIndex() const = 0;

   /** As we iterate along, we keep a column map entry for each sequence
    * visited.  This works out pretty well except for extreme cases (such
    * as iterating over entire fly genomes where we can accumulate 10s of 
    * thousands of empty entries for all the different scaffolds when 
    * in truth we only need a handful at any given time). Under these
    * circumstances, calling this method every 1M bases or so will help
    * reduce memory as well as speed up queries on the column map. Perhaps
    * this should eventually be built in and made transparent? */
   virtual void defragment() const = 0;
   
protected:
   friend class counted_ptr<ColumnIterator>;
   friend class counted_ptr<const ColumnIterator>;
   virtual ~ColumnIterator() = 0;
};

inline ColumnIterator::~ColumnIterator() {}

inline bool operator<(ColumnIteratorConstPtr p1,
                      ColumnIteratorConstPtr p2) 
{
  if (p1.get() == NULL && p2.get() == NULL)
  {
    return false;
  }
  else if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL;
  }
  else
  {
    return p1->leftOf(p2);
  }
}

}

#endif
