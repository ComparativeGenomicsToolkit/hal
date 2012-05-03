/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTCOLUMNITERATOR_H
#define _DEFAULTCOLUMNITERATOR_H

#include <set>
#include "halColumnIterator.h"

namespace hal {

class DefaultColumnIterator : public ColumnIterator
{
public:

   DefaultColumnIterator(hal::Genome* reference, hal::Genome* root = NULL,
                         hal_index_t columnIndex = 0);
   
   ~DefaultColumnIterator();

   /** Move column iterator one column to the right along reference
    * genoem sequence */
    void toRight() const;

   /** Test if column iterator is at the same position (only in terms
    * of the reference genome) as another iterator.  Used to bound
    * a loop, for example */
    bool equals(ColumnIteratorConstPtr other) const;
   
   /** Direct read access to the segment iterator map (for lack of a 
    * better interface for now) */
   const SegmentMap& getSegmentMap() const;

   const hal::Genome* getReferenceGenome() const;

private:
   
   typedef std::pair<const Genome*, hal_index_t> VisitFlag;
   typedef std::set<VisitFlag> VisitSet;
   
   VisitSet _topVisited;
   VisitSet _bottomVisited;

   SegmentMap _segmentMap;
   
};

}

#endif
