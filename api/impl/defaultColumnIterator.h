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
#include "halRearrangement.h"

namespace hal {

class DefaultColumnIterator : public ColumnIterator
{
public:

   DefaultColumnIterator(const hal::Sequence* reference, 
                         const hal::Genome* root,
                         hal_index_t columnIndex,
                         hal_index_t lastIndex,
                         hal_size_t maxInsertionLength,
                         bool noDupes);
   
   ~DefaultColumnIterator();

   /** Move column iterator one column to the right along reference
    * genoem sequence */
    void toRight() const;

   bool lastColumn() const;

   const hal::Genome* getReferenceGenome() const;
   const hal::Sequence* getReferenceSequence() const;

   /** Get a pointer to the column map */
   const ColumnMap* getColumnMap() const;

   hal_index_t getArrayIndex() const;

   void defragment() const;

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
   typedef std::map<const Genome*, VisitSet*> VisitCache;

   struct StackEntry 
   {
      const Sequence* _sequence;
      hal_index_t _firstIndex;
      hal_index_t _index;
      hal_index_t _lastIndex;
      hal_size_t _cumSize; 
      LinkedTopIteratorPtr _top;
      LinkedBottomIteratorPtr _bottom;
      RearrangementPtr _rearrangement;
   };

   typedef std::vector<StackEntry> ActiveStack;

private:
   void pushStack(ActiveStack& stack, const Sequence* ref, hal_index_t index, 
                  hal_index_t lastIndex, bool update) const;
   bool handleDeletion(TopSegmentIteratorConstPtr inputTopIterator) const;
   bool handleInsertion(TopSegmentIteratorConstPtr inputTopIterator) const;
   void resetColMap() const;
   void eraseColMap() const;
   void recursiveUpdate(bool init) const;

   void updateParent(LinkedTopIteratorPtr topIt) const;
   void updateChild(LinkedBottomIteratorPtr bottomIt, 
                    hal_size_t index) const;
   void updateNextTopDup(LinkedTopIteratorPtr topIt) const;
   void updateParseUp(LinkedBottomIteratorPtr bottomIt) const;
   void updateParseDown(LinkedTopIteratorPtr topIt) const;

   bool inBounds() const;
   bool nextFreeIndex() const;
   bool colMapInsert(DNAIteratorConstPtr dnaIt) const;
   bool checkRange(DNAIteratorConstPtr dnaIt) const;
   
private:

   // everything's mutable to keep const behaviour consistent with
   // other iterators (which provide both const and non-const access)
   // the fact that this iterator has no writable interface makes it
   // seem like a dumb excercise though. 
   mutable const Genome* _root;
   mutable ActiveStack _stack;
   mutable ActiveStack _indelStack;
   mutable const Sequence* _ref;
   mutable size_t _curInsertionLength;

   mutable hal_size_t _maxInsertionLength;
   mutable bool _noDupes;

   mutable ColumnMap _colMap;
   mutable TopSegmentIteratorConstPtr _top;
   mutable TopSegmentIteratorConstPtr _next;
   mutable VisitCache _visitCache;
};




}
#endif

