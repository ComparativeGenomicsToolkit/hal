/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNITERATOR_H
#define _HALCOLUMNITERATOR_H

#include "halColumnIteratorStack.h"
#include "halDefs.h"
#include "halDnaIterator.h"
#include "halPositionCache.h"
#include "halSequence.h"
#include "sonLib.h"
#include <list>
#include <map>
#include <set>

namespace hal {

    /**
     * Interface Column iterator for allowing traditional maf-like (left-to-right)
     * parsing of a hal alignment.  Columns are iterated with respect to
     * a specified reference genome.  This isn't the most efficient way
     * to explore the hal structure, which is designed for bottom-up and/or
     * top-down traversal.
     */
    class ColumnIterator {
      public:
        /* constructor */
        ColumnIterator(const Genome *reference, const std::set<const Genome *> *targets, hal_index_t columnIndex,
                       hal_index_t lastColumnIndex, hal_size_t maxInsertLength, bool noDupes, bool noAncestors,
                       bool reverseStrand, bool unique, bool onlyOrthologs);

        /** Destructor */
        virtual ~ColumnIterator();

        /// @cond TEST
        // Originally could compared genomes by pointers (because they are
        // persistent and unique).  However this lead to output instability
        // problems for tests, so we use names.  If this shows up as a performance issue,
        // we can sort only in MAF code.
        struct SequenceLess {
            bool operator()(const Sequence *s1, const Sequence *s2) const {
                int diff = s1->getGenome()->getName().compare(s2->getGenome()->getName());
                return (diff < 0) || ((diff == 0) && s1->getArrayIndex() < s2->getArrayIndex());
            }
        };
        /// @endcond

        typedef std::vector<DnaIteratorPtr> DNASet;
        typedef std::map<const Sequence *, DNASet *, SequenceLess> ColumnMap;

        /** Move column iterator one column to the right along reference
         * genoem sequence */
        virtual void toRight();

        /** Move column iterator to arbitrary site in genome -- effectively
         * resetting the iterator (convenience function to avoid creation of
         * new iterators in some cases).
         * @param columnIndex position of column in forward genome coordinates
         * @param lastIndex last column position (for iteration).  must be greater
         *  than columnIndex
         * @param clearCache clear the cache that prevents columns from being
         * visited twice.  If not set to true, then its possible the iterator
         * ends up not at "columnIndex" but at the next unvisited column.*/
        virtual void toSite(hal_index_t columnIndex, hal_index_t lastIndex, bool clearCache = false);

        /** Use this method to bound iteration loops.  When the column iterator
         * is retrieved from the sequence or genome, the last column is specified.
         * toRight() can then be called until lastColumn is true.  */
        virtual bool lastColumn() const;

        /** Get a pointer to the reference genome for the column iterator */
        virtual const Genome *getReferenceGenome() const;

        /** Get a pointer to the reference sequence for the column iterator */
        virtual const Sequence *getReferenceSequence() const;

        /** Get the position in the reference sequence
         * NOTE
         * Seems to be returning the next position, rather than the current.
         * Must go back and review but it is concerning. */
        virtual hal_index_t getReferenceSequencePosition() const;

        /** Get a pointer to the column map */
        virtual const ColumnMap *getColumnMap() const;

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
        virtual void defragment();

        /** Check whether the column iterator's left-most reference coordinate
         * is within the iterator's range, ie is "canonical".  This can be used
         * to ensure that the same reference position does not get sampled by
         * different iterators covering distinct ranges.  If there are no
         * duplications, then this function will always return true. */
        virtual bool isCanonicalOnRef() const;

        /** Print contents of column iterator */
        virtual void print(std::ostream &os) const;

        /** Get a new tree that represents the phylogenetic relationship
         * between the entries in this column. Do not attempt to free this
         * tree. */
        virtual stTree *getTree() const;

        // temp -- probably want to have a "global column iterator" object
        // instead
        typedef std::map<const Genome *, PositionCache *> VisitCache;
        virtual VisitCache *getVisitCache();
        virtual void setVisitCache(VisitCache *visitCache);
        virtual void clearVisitCache();

      private:
        typedef ColumnIteratorStack::LinkedBottomIterator LinkedBottomIterator;
        typedef ColumnIteratorStack::LinkedTopIterator LinkedTopIterator;
        typedef ColumnIteratorStack::Entry StackEntry;

      private:
        void recursiveUpdate(bool init);
        bool handleDeletion(const TopSegmentIteratorPtr &inputTopSegIt);
        bool handleInsertion(const TopSegmentIteratorPtr &inputTopSegIt);

        void updateParent(LinkedTopIterator *linkTopIt);
        void updateChild(LinkedBottomIterator *linkBotIt, hal_size_t index);
        void updateNextTopDup(LinkedTopIterator *linkTopIt);
        void updateParseUp(LinkedBottomIterator *linkBotIt);
        void updateParseDown(LinkedTopIterator *linkTopIt);

        bool parentInScope(const Genome *) const;
        bool childInScope(const Genome *, hal_size_t child) const;
        void nextFreeIndex();
        bool colMapInsert(DnaIteratorPtr dnaIt);

        void resetColMap();
        void eraseColMap();

        stTree *buildTree() const;
        void clearTree();

      private:
        std::set<const Genome *> _targets;
        std::set<const Genome *> _scope;
        ColumnIteratorStack _stack;
        ColumnIteratorStack _indelStack;
        const Sequence *_refSeq;

        RearrangementPtr _rearrangement;
        hal_size_t _maxInsertionLength;
        bool _noDupes;
        bool _noAncestors;

        ColumnMap _colMap;
        TopSegmentIteratorPtr _top;
        TopSegmentIteratorPtr _next;
        VisitCache _visitCache;
        bool _break;
        const Sequence *_prevRefSequence;
        hal_index_t _prevRefIndex;
        hal_index_t _leftmostRefPos;
        mutable stTree *_treeCache;
        bool _unique;
        bool _onlyOrthologs;
    };

    inline std::ostream &operator<<(std::ostream &os, const ColumnIterator &cit) {
        cit.print(os);
        return os;
    }
    inline bool ColumnIterator::parentInScope(const Genome *genome) const {
        assert(genome != NULL && genome->getParent() != NULL);
        return _scope.empty() || _scope.find(genome->getParent()) != _scope.end();
    }

    inline bool ColumnIterator::childInScope(const Genome *genome, hal_size_t child) const {
        assert(genome != NULL && genome->getChild(child) != NULL);
        return _scope.empty() || _scope.find(genome->getChild(child)) != _scope.end();
    }
}

#endif
// Local Variables:
// mode: c++
// End:
