/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNITERATOR_H
#define _HALCOLUMNITERATOR_H

#include <set>
#include <list>
#include <map>
#include <set>
#include "sonLib.h"
#include "halPositionCache.h"
#include "halDefs.h"
#include "halDNAIterator.h"
#include "halSequence.h"
#include "halColumnIteratorStack.h"

namespace hal {

/** 
 * Interface Column iterator for allowing traditional maf-like (left-to-right)  
 * parsing of a hal alignment.  Columns are iterated with respect to
 * a specified reference genome.  This isn't the most efficient way
 * to explore the hal structure, which is designed for bottom-up and/or
 * top-down traversal.  
 */
class ColumnIterator 
{
public:

    /* constructor */
    ColumnIterator(const Genome* reference, 
                   const std::set<const Genome*>* targets,
                   hal_index_t columnIndex,
                   hal_index_t lastColumnIndex,
                   hal_size_t maxInsertLength,
                   bool noDupes,
                   bool noAncestors,
                   bool reverseStrand,
                   bool unique,
                   bool onlyOrthologs);
    
   /// @cond TEST
   // we can compare genomes by pointers (because they are persistent
   // and unique, though it's still hacky) but we can't do the same 
   // for anything else, including sequences.  
   struct SequenceLess { bool operator()(const Sequence* s1,
                                         const Sequence* s2) const {
     return s1->getGenome() < s2->getGenome() || (
       s1->getGenome() == s2->getGenome() && 
       s1->getArrayIndex() < s2->getArrayIndex()); }
   };
   /// @endcond

   typedef std::vector<DNAIteratorConstPtr> DNASet;
   typedef std::map<const Sequence*, DNASet*, SequenceLess> ColumnMap;

   /** Move column iterator one column to the right along reference
    * genoem sequence */
   virtual void toRight() const;

   /** Move column iterator to arbitrary site in genome -- effectively
    * resetting the iterator (convenience function to avoid creation of
    * new iterators in some cases).  
    * @param columnIndex position of column in forward genome coordinates 
    * @param lastIndex last column position (for iteration).  must be greater
    *  than columnIndex 
    * @param clearCache clear the cache that prevents columns from being 
    * visited twice.  If not set to true, then its possible the iterator
    * ends up not at "columnIndex" but at the next unvisited column.*/
   virtual void toSite(hal_index_t columnIndex, 
                       hal_index_t lastIndex,
                       bool clearCache = false) const;

   /** Use this method to bound iteration loops.  When the column iterator
    * is retrieved from the sequence or genome, the last column is specfied.
    * toRight() cna then be called until lastColumn is true.  */
   virtual bool lastColumn() const;
   
   /** Get a pointer to the reference genome for the column iterator */
   virtual const Genome* getReferenceGenome() const;

   /** Get a pointer to the reference sequence for the column iterator */
   virtual const Sequence* getReferenceSequence() const;

   /** Get the position in the reference sequence 
    * NOTE
    * Seems to be returning the next position, rather than the current.
    * Must go back and review but it is concerning. */
   virtual hal_index_t getReferenceSequencePosition() const;

   /** Get a pointer to the column map */
   virtual const ColumnMap* getColumnMap() const;

   /** Get the index of the column in the reference genome's array */
   virtual hal_index_t getArrayIndex() const;

   /** As we iterate along, we keep a column map entry for each sequence
    * visited.  This works out pretty well except for extreme cases (such
    * as iterating over entire fly genomes where we can accumulate 10s of 
    * thousands of empty entries for all the different scaffolds when 
    * in truth we only need a handful at any given time). Under these
    * circumstances, calling this method every 1M bases or so will help
    * reduce memory as well as speed up queries on the column map. Perhaps
    * this should eventually be built in and made transparent? */
   virtual void defragment() const;

   /** Check whether the column iterator's left-most reference coordinate
    * is within the iterator's range, ie is "canonical".  This can be used
    * to ensure that the same reference position does not get sampled by
    * different iterators covering distinct ranges.  If there are no 
    * duplications, then this function will always return true. */
   virtual bool isCanonicalOnRef() const;

   /** Print contents of column iterator */
   virtual void print(std::ostream& os) const;

   /** Get a new tree that represents the phylogenetic relationship
    * between the entries in this column. Do not attempt to free this
    * tree. */
   virtual stTree *getTree() const;

   // temp -- probably want to have a "global column iterator" object
   // instead
   typedef std::map<const Genome*, PositionCache*> VisitCache;
   virtual VisitCache *getVisitCache() const;
   virtual void setVisitCache(VisitCache *visitCache) const;
   virtual void clearVisitCache() const;

protected:

   typedef ColumnIteratorStack::LinkedBottomIterator LinkedBottomIterator;
   typedef ColumnIteratorStack::LinkedTopIterator LinkedTopIterator;
   typedef ColumnIteratorStack::Entry StackEntry;
private:

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

private:

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

private:
   friend class counted_ptr<ColumnIterator>;
   friend class counted_ptr<const ColumnIterator>;
   virtual ~ColumnIterator();
};

inline std::ostream& operator<<(std::ostream& os, const ColumnIterator& cit)
{
  cit.print(os);
  return os;
}
inline bool ColumnIterator::parentInScope(const Genome* genome) const
{
  assert(genome != NULL && genome->getParent() != NULL);
  return _scope.empty() || _scope.find(genome->getParent()) != _scope.end();
}

inline bool ColumnIterator::childInScope(const Genome* genome,
                                         hal_size_t child) const
{
  assert(genome != NULL && genome->getChild(child) != NULL);
  return _scope.empty() || _scope.find(genome->getChild(child)) != _scope.end();
}

}


#endif
