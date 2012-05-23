/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTCOLUMNITERATOR_H
#define _DEFAULTCOLUMNITERATOR_H

#include <set>
#include <stack>
#include <vector>
#include <map>
#include "halColumnIterator.h"

namespace hal {

class DefaultColumnIterator : public ColumnIterator
{
public:

   DefaultColumnIterator(const hal::Genome* reference, 
                         const hal::Genome* root = NULL,
                         hal_index_t columnIndex = 0,
                         hal_size_t maxInsertionLength = 500);
   
   ~DefaultColumnIterator();

   /** Move column iterator one column to the right along reference
    * genoem sequence */
    void toRight() const;

    bool leftOf(ColumnIteratorConstPtr other) const;
   
   const hal::Genome* getReferenceGenome() const;

   /** Get a pointer to the column map */
   const ColumnMap* getColumnMap() const;

private:

   struct LinkedTopIterator;
   typedef smart_ptr<LinkedTopIterator>::type LinkedTopIteratorPtr;
   struct LinkedBottomIterator;
   typedef smart_ptr<LinkedBottomIterator>::type LinkedBottomIteratorPtr;

   struct LinkedBottomIterator 
   {
      BottomSegmentIteratorConstPtr _it;
      DNAIteratorConstPtr _dna;
      LinkedTopIteratorPtr _topParse;
      std::vector<LinkedTopIteratorPtr> _children;
      LinkedBottomIteratorPtr _nextDup;
   };

   struct LinkedTopIterator 
   {
      TopSegmentIteratorConstPtr _it;
      DNAIteratorConstPtr _dna;
      LinkedBottomIteratorPtr _bottomParse;
      LinkedBottomIteratorPtr _parent;
      LinkedTopIteratorPtr _nextDup;
   };

   typedef std::pair<const Genome*, hal_index_t> VisitFlag;
   typedef std::set<VisitFlag> VisitSet;
   typedef std::map<const Genome*, LinkedBottomIteratorPtr> BottomMap;
   typedef std::map<const Genome*, LinkedTopIteratorPtr> TopMap;
   typedef std::stack<const Genome*> GenomeStack;

private:

   void init() const;
   void resetColMap() const;
   hal_index_t moveRightToNextUnvisited(LinkedTopIteratorPtr topIt) const;

   void updateParent(LinkedTopIteratorPtr topIt) const;
   void updateChild(LinkedBottomIteratorPtr bottomIt, 
                    hal_size_t index) const;
   void updateNextTopDup(LinkedTopIteratorPtr topIt) const;
   void updateNextBottomDup(LinkedBottomIteratorPtr bottomIt) const;
   void updateParseUp(LinkedBottomIteratorPtr bottomIt) const;
   void updateParseDown(LinkedTopIteratorPtr topIt) const;
   
private:

   // everything's mutable to keep const behaviour consistent with
   // other iterators (which provide both const and non-const access)
   // the fact that this iterator has no writable interface makes it
   // seem like a dumb excercise though. 
   mutable VisitSet _topVisited;
   mutable VisitSet _bottomVisited;

   mutable const Genome* _root;
   mutable const Genome* _reference;
   mutable std::stack<const Genome*> _activeStack;
   mutable size_t _curInsertionLength;

   mutable BottomMap _bottomMap;
   mutable TopMap _topMap;

   mutable hal_index_t _index;
   mutable hal_size_t _maxInsertionLength;

   mutable ColumnMap _colMap;
};

}

#endif

