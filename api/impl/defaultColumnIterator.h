/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTCOLUMNITERATOR_H
#define _DEFAULTCOLUMNITERATOR_H

#include <set>
#include <stack>
#include "halColumnIterator.h"

#ifdef _DISABLE_FOR_NOW_
namespace hal {

class DefaultColumnIterator : public ColumnIterator
{
public:

   DefaultColumnIterator(hal::Genome* reference, hal::Genome* root = NULL,
                         hal_index_t columnIndex = 0,
                         hal_size_t maxInsertionLength = 500);
   
   ~DefaultColumnIterator();

   /** Move column iterator one column to the right along reference
    * genoem sequence */
    void toRight() const;

   /** Test if column iterator is at the same position (only in terms
    * of the reference genome) as another iterator.  Used to bound
    * a loop, for example */
    bool equals(ColumnIteratorConstPtr other) const;
   
   const hal::Genome* getReferenceGenome() const;

private:

   void init() const;
   void update() const;

private:
   
   typedef std::pair<const Genome*, hal_index_t> VisitFlag;
   typedef std::set<VisitFlag> VisitSet;
   typedef std::map<const Genome*, LinkedBottomIterator> BottomMap;
   typedef std::map<const Genome*, LinkedTopIterator> TopMap;

   struct LinkedTopIterator;
   struct LinkedBottomIterator 
   {
      BottomIteratorConstPtr _it;
      LinkedTopIterator _topParse;
      std::vector<LinkedTopIterator> _children;
      LinkedBottomIterator _nextDup;
   };

   struct LinkedTopIterator 
   {
      TopIteratorConstPtr _it;
      LinkedBottomIterator _bottomParse;
      LinkedBottomIterator _parent;
      LinkedTopIterator _nextDup;
   };

   // everything's mutable to keep const behaviour consistent with
   // other iterators (which provide both const and non-const access)
   // the fact that this iterator has no writable interface makes it
   // seem like a dumb excercise though. 
   mutable VisitSet _topVisited;
   mutable VisitSet _bottomVisited;

   mutable const Genome* _root;
   mutable std::stack<const Genome*> _activeStack;

   mutable BottomMap _bottomMap;
   mutable TopMap _topMap;

   mutable hal_index_t _index;
   mutable hal_size_t _maxInsertionLength;
};

}

#endif
#endif
