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
#include "halCommon.h"
#include "columnIteratorStack.h"

namespace hal {

class DefaultColumnIterator : public ColumnIterator
{
public:
   DefaultColumnIterator(const hal::Genome* reference, 
                         const std::set<const Genome*>* targets,
                         hal_index_t columnIndex,
                         hal_index_t lastIndex,
                         hal_size_t maxInsertionLength,
                         bool noDupes,
                         bool noAncestors,
                         bool reverseStrand,
                         bool unique,
                         bool onlyOrthologs);
   
   virtual ~DefaultColumnIterator();

   // COLUMN ITERATOR INTERFACE
   virtual void toRight() const;
   virtual void toSite(hal_index_t columnIndex, hal_index_t lastIndex,
                       bool clearCache) const;
   virtual bool lastColumn() const;
   virtual const hal::Genome* getReferenceGenome() const;
   virtual const hal::Sequence* getReferenceSequence() const;
   virtual hal_index_t getReferenceSequencePosition() const;
   virtual const ColumnMap* getColumnMap() const;
   virtual hal_index_t getArrayIndex() const;
   virtual void defragment() const;
   virtual bool isCanonicalOnRef() const;
   virtual void print(std::ostream& os) const;
   virtual stTree *getTree() const;
   virtual VisitCache *getVisitCache() const;
   virtual void setVisitCache(VisitCache *visitCache) const;
   virtual void clearVisitCache() const;
protected:

   typedef ColumnIteratorStack::LinkedBottomIterator LinkedBottomIterator;
   typedef ColumnIteratorStack::LinkedTopIterator LinkedTopIterator;
   typedef ColumnIteratorStack::Entry StackEntry;
   
protected:

   void recursiveUpdate(bool init) const;
   bool handleDeletion(TopSegmentIteratorConstPtr inputTopIterator) const;
   bool handleInsertion(TopSegmentIteratorConstPtr inputTopIterator) const;

   void updateParent(LinkedTopIterator* topIt) const;
   void updateChild(LinkedBottomIterator* bottomIt, hal_size_t index) const;
   void updateNextTopDup(LinkedTopIterator* topIt) const;
   void updateParseUp(LinkedBottomIterator* bottomIt) const;
   void updateParseDown(LinkedTopIterator* topIt) const;

   bool parentInScope(const Genome*) const;
   bool childInScope(const Genome*, hal_size_t child) const;
   void nextFreeIndex() const;
   bool colMapInsert(DNAIteratorConstPtr dnaIt) const;

   void resetColMap() const;
   void eraseColMap() const;

   void clearTree() const;

protected:

   // everything's mutable to keep const behaviour consistent with
   // other iterators (which provide both const and non-const access)
   // the fact that this iterator has no writable interface makes it
   // seem like a dumb excercise though. 
   mutable std::set<const Genome*> _targets;
   mutable std::set<const Genome*> _scope;
   mutable ColumnIteratorStack _stack;
   mutable ColumnIteratorStack _indelStack;
   mutable const Sequence* _ref;
   mutable size_t _curInsertionLength;

   mutable RearrangementPtr _rearrangement;
   mutable hal_size_t _maxInsertionLength;
   mutable bool _noDupes;
   mutable bool _noAncestors;
   mutable bool _reversed;

   mutable ColumnMap _colMap;
   mutable TopSegmentIteratorConstPtr _top;
   mutable TopSegmentIteratorConstPtr _next;
   mutable VisitCache _visitCache;
   mutable bool _break;
   mutable const Sequence* _prevRefSequence;
   mutable hal_index_t _prevRefIndex;
   mutable hal_index_t _leftmostRefPos;
   mutable stTree *_tree;
   mutable bool _unique;
   mutable bool _onlyOrthologs;
};

inline bool DefaultColumnIterator::parentInScope(const Genome* genome) const
{
  assert(genome != NULL && genome->getParent() != NULL);
  return _scope.empty() || _scope.find(genome->getParent()) != _scope.end();
}

inline bool DefaultColumnIterator::childInScope(const Genome* genome,
                                                hal_size_t child) const
{
  assert(genome != NULL && genome->getChild(child) != NULL);
  return _scope.empty() || _scope.find(genome->getChild(child)) != _scope.end();
}

}
#endif

