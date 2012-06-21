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

   typedef std::set<hal::DNAIteratorConstPtr> DNASet;
   typedef std::map<const hal::Sequence*, DNASet> ColumnMap;

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
