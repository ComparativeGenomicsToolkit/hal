/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _COLUMNITERATORSTACK_H
#define _COLUMNITERATORSTACK_H

#include "halCommon.h"
#include "halRearrangement.h"
#include <map>
#include <set>
#include <stack>
#include <vector>

namespace hal {

    class ColumnIteratorStack {
      public:
        class Entry;
        class LinkedTopIterator;
        class LinkedBottomIterator {
          public:
            LinkedBottomIterator() : _topParse(NULL), _entry(NULL) {
            }
            BottomSegmentIteratorPtr _it;
            DnaIteratorPtr _dna;
            LinkedTopIterator *_topParse;
            std::vector<LinkedTopIterator *> _children;
            Entry *_entry;
        };

        class LinkedTopIterator {
          public:
            LinkedTopIterator() : _bottomParse(NULL), _parent(NULL), _nextDup(NULL), _entry(NULL) {
            }
            TopSegmentIteratorPtr _it;
            DnaIteratorPtr _dna;
            LinkedBottomIterator *_bottomParse;
            LinkedBottomIterator *_parent;
            LinkedTopIterator *_nextDup;
            Entry *_entry;
        };

        class Entry {
          public:
            Entry(const Sequence *seq, hal_index_t first, hal_index_t index, hal_index_t last, hal_size_t size, bool reversed)
                : _sequence(seq), _firstIndex(first), _index(index), _lastIndex(last), _cumulativeSize(size), _reversed(reversed) {
                _top._entry = this;
                _bottom._entry = this;
            }
            ~Entry() {
                freeLinks();
            }

            bool pastEnd() const {
                if (_reversed) {
                    return _index < _firstIndex;
                } else {
                    return _index > _lastIndex;
                }
            }

            LinkedTopIterator *newTop() {
                LinkedTopIterator *top = new LinkedTopIterator();
                top->_entry = this;
                _topLinks.push_back(top);
                return top;
            }

            LinkedBottomIterator *newBottom() {
                LinkedBottomIterator *bottom = new LinkedBottomIterator();
                bottom->_entry = this;
                _bottomLinks.push_back(bottom);
                return bottom;
            }

            void freeLinks() {
                size_t i;
                for (i = 0; i < _topLinks.size(); ++i) {
                    delete _topLinks[i];
                }
                _topLinks.clear();
                _top._bottomParse = NULL;
                _top._parent = NULL;
                _top._nextDup = NULL;

                for (i = 0; i < _bottomLinks.size(); ++i) {
                    delete _bottomLinks[i];
                }
                _bottomLinks.clear();
                _bottom._topParse = NULL;
                _bottom._children.clear();
            }
            const Sequence *_sequence;
            hal_index_t _firstIndex;
            hal_index_t _index;
            hal_index_t _lastIndex;
            hal_size_t _cumulativeSize;
            LinkedTopIterator _top;
            LinkedBottomIterator _bottom;
            std::vector<LinkedTopIterator *> _topLinks;
            std::vector<LinkedBottomIterator *> _bottomLinks;
            bool _reversed;
        };

      public:
        ~ColumnIteratorStack() {
            clear();
        }
        void push(const Sequence *ref, hal_index_t index, hal_index_t lastIndex, bool reversed = false) {
            assert(lastIndex >= index);
            assert(ref != NULL);
            hal_size_t cumulative = 0;
            if (_stack.size() > 0) {
                cumulative = top()->_cumulativeSize + lastIndex - index + 1;
            }
            Entry *entry = new Entry(ref, index, reversed ? lastIndex : index, lastIndex, cumulative, reversed);
            _stack.push_back(entry);
        }
        void pushStack(ColumnIteratorStack &otherStack) {
            for (size_t i = 0; i < otherStack.size(); ++i) {
                _stack.push_back(otherStack[i]);
            }
            otherStack._stack.clear();
        }
        void popDelete() {
            delete _stack.back();
            _stack.pop_back();
        }
        void clear() {
            while (_stack.size() > 0) {
                popDelete();
            }
        }
        const Entry *top() const {
            return _stack.back();
        }
        Entry *top() {
            return _stack.back();
        }
        const Entry *operator[](size_t i) const {
            return _stack[i];
        }
        Entry *operator[](size_t i) {
            return _stack[i];
        }
        size_t size() const {
            return _stack.size();
        }
        void resetLinks() {
            for (size_t i = 0; i < _stack.size(); ++i) {
                _stack[i]->freeLinks();
            }
        }
        bool topInBounds() const {
            const Entry *e = _stack.back();
            return e->_index >= e->_firstIndex && e->_index <= e->_lastIndex;
        }

      private:
        std::vector<Entry *> _stack;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
