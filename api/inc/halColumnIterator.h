/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNITERATOR_H
#define _HALCOLUMNITERATOR_H

#include <list>
#include <map>
#include "halDefs.h"

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

   typedef std::vector<hal::DNAIteratorConstPtr> DNAList;
   typedef std::map<const hal::Genome*, DNAList> ColumnMap;

   /** Move column iterator one column to the right along reference
    * genoem sequence */
   virtual void toRight() const = 0;

   /** Test if column iterator is at the same position (only in terms
    * of the reference genome) as another iterator.  Used to bound
    * a loop, for example */
   virtual bool equals(ColumnIteratorConstPtr other) const = 0;
   
   /** Get a pointer to the reference genome for the column iterator */
   virtual const hal::Genome* getReferenceGenome() const = 0;

   /** Get a pointer to the column map */
   virtual const ColumnMap* getColumnMap() const = 0;
   
protected:
   friend class counted_ptr<ColumnIterator>;
   friend class counted_ptr<const ColumnIterator>;
   virtual ~ColumnIterator() = 0;
};

inline ColumnIterator::~ColumnIterator() {}

inline bool operator==(ColumnIteratorConstPtr p1,
                       ColumnIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(ColumnIteratorConstPtr p1,
                       ColumnIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}

#endif
