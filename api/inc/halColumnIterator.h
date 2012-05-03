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

   /** Linked list of top segment iterators */
   typedef std::list<TopSegmentIteratorConstPtr> TopList;
   /** Linked list of bottom segment iterators */
   typedef std::list<BottomSegmentIteratorConstPtr> BottomList;
   /** Pair of iterators */
   typedef std::pair<TopList, BottomList> ListPair;
   /** Map (bijective) between set of genomes and lists of iterators 
    * the first iterator (pair) is the current segment for the genome.
    * all remaining pairs of duplications of that segment */
   typedef std::map<const Genome*, ListPair> SegmentMap;

   /** Move column iterator one column to the right along reference
    * genoem sequence */
   virtual void toRight() const = 0;

   /** Test if column iterator is at the same position (only in terms
    * of the reference genome) as another iterator.  Used to bound
    * a loop, for example */
   virtual bool equals(ColumnIteratorConstPtr other) const = 0;
   
   /** Direct read access to the segment iterator map (for lack of a 
    * better interface for now) */
   virtual const SegmentMap& getSegmentMap() = 0;

   /** Get a pointer to the reference genome for the column iterator */
   virtual const hal::Genome* getReferenceGenome() const = 0;
   
protected:
   friend class counted_ptr<ColumnIterator>;
   friend class counted_ptr<const ColumnIterator>;
   virtual ~ColumnIterator() = 0;
};

}

#endif
