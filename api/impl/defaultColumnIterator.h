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

   DefaultColumnIterator(const hal::Sequence* reference, 
                         const hal::Genome* root,
                         hal_index_t columnIndex,
                         hal_size_t maxInsertionLength,
                         bool endIterator);
   
   ~DefaultColumnIterator();

   /** Move column iterator one column to the right along reference
    * genoem sequence */
    void toRight() const;

    bool leftOf(ColumnIteratorConstPtr other) const;
   
   const hal::Genome* getReferenceGenome() const;
   const hal::Sequence* getReferenceSequence() const;

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

   typedef std::set<hal_index_t> VisitSet;

   struct StackEntry 
   {
      const Sequence* _sequence;
      hal_index_t _index;
      bool _reversed;
      LinkedTopIteratorPtr _top;
      LinkedBottomIteratorPtr _bottom;
      VisitSet _visitSet;
   };

   typedef std::stack<StackEntry> ActiveStack;

private:

   void init(const hal::Sequence* ref, hal_index_t index, 
             bool endIterator) const;
   void resetColMap() const;
   void eraseColMap() const;
   void recursiveUpdate(bool init) const;
   hal_index_t moveRightToNextUnvisited(LinkedTopIteratorPtr topIt) const;

   void updateParent(LinkedTopIteratorPtr topIt) const;
   void updateChild(LinkedBottomIteratorPtr bottomIt, 
                    hal_size_t index) const;
   void updateNextTopDup(LinkedTopIteratorPtr topIt) const;
   void updateNextBottomDup(LinkedBottomIteratorPtr bottomIt) const;
   void updateParseUp(LinkedBottomIteratorPtr bottomIt) const;
   void updateParseDown(LinkedTopIteratorPtr topIt) const;

   hal_index_t nextFreeIndex(LinkedTopIteratorPtr topIt) const;
   void colMapInsert(DNAIteratorConstPtr dnaIt, 
                     bool updateVisitSet = true) const;
   bool checkRange(DNAIteratorConstPtr dnaIt) const;

   
private:

   // everything's mutable to keep const behaviour consistent with
   // other iterators (which provide both const and non-const access)
   // the fact that this iterator has no writable interface makes it
   // seem like a dumb excercise though. 
   mutable const Genome* _root;
   mutable ActiveStack _stack;
   mutable size_t _curInsertionLength;

   mutable hal_size_t _maxInsertionLength;

   mutable ColumnMap _colMap;
};




}
#endif

