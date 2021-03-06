/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halRearrangement.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedTopSegmentIterator.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

// maximum size a simple indel can be to be considered a gap (and not
// a rearrangement)
const hal_size_t Rearrangement::DefaultGapThreshold = 10;
const double Rearrangement::DefaultNThreshold = 0.10;

Rearrangement::Rearrangement(const Genome *childGenome, hal_size_t gapThreshold, double nThreshold, bool atomic)
    : _gapThreshold(gapThreshold), _atomic(atomic), _nThreshold(nThreshold), _genome(childGenome), _parent(NULL), _id(Invalid) {
    _parent = childGenome->getParent();
    assert(_parent != NULL);
    // just allocating here -- need to be properly init
    _cur = _genome->getGappedTopSegmentIterator(0, _gapThreshold, _atomic);
    assert(_cur->getGapThreshold() == _gapThreshold);
    _next = _cur->clone();
    _left = _cur->clone();
    _right = _left->clone();
    _leftParent = _parent->getGappedBottomSegmentIterator(0, 0, _gapThreshold, _atomic);
    _rightParent = _leftParent->clone();
    _curParent = _leftParent->clone();
    _nextParent = _leftParent->clone();
    _top = _cur->getLeft()->clone();
}

// Rearrangement Interface Methods
Rearrangement::ID Rearrangement::getID() const {
    return _id;
}

hal_size_t Rearrangement::getLength() const {
    return _id == Deletion ? _leftParent->getLength() : _cur->getLength();
}

hal_size_t Rearrangement::getNumContainedGaps() const {
    return _id == Deletion ? _leftParent->getNumGaps() : _cur->getNumGaps();
}

hal_size_t Rearrangement::getNumContainedGapBases() const {
    return _id == Deletion ? _leftParent->getNumGapBases() : _cur->getNumGapBases();
}

TopSegmentIteratorPtr Rearrangement::getLeftBreakpoint() const {
    return _cur->getLeft();
}

TopSegmentIteratorPtr Rearrangement::getRightBreakpoint() const {
    return _cur->getRight();
}

bool Rearrangement::identifyFromLeftBreakpoint(TopSegmentIteratorPtr topSegment) {
    if (scanDuplicationCycle(topSegment) == true) {
        _id = _cur->isMissingData(_nThreshold) ? Missing : Duplication;
    } else if (scanNothingCycle(topSegment) == true) {
        _id = Nothing;
    } else if (scanInversionCycle(topSegment) == true) {
        _id = _cur->isMissingData(_nThreshold) ? Missing : Inversion;
    } else if (scanInsertionCycle(topSegment) == true) {
        if (_cur->isMissingData(_nThreshold) == true) {
            _id = Missing;
        } else if (_cur->hasParent() == true && !_cur->isFirst() && !_cur->isLast()) {
            _id = Transposition;
        } else if (_cur->hasParent() == false) {
            _id = _cur->getLength() > _gapThreshold ? Insertion : Gap;
        }
    } else if (scanDeletionCycle(topSegment) == true) {
        if (_leftParent->isMissingData(_nThreshold) == true) {
            _id = Missing;
        } else if (_leftParent->hasChild() == false) {
            _id = _leftParent->getLength() > _gapThreshold ? Deletion : Gap;
        } else {
            _id = Complex;
        }
    } else {
        /*
            resetStatus(topSegment);
            if (_cur->isFirst() == false && _cur->isLast() == false)
            {
              _left->toLeft();
              _right->toRight();
              if (_left->hasParent() == true && _right->hasParent() == true &&
                  _cur->hasParent() == true)
              {
                _curParent->toParent(_cur);
                _leftParent->toParent(_left);
                _rightParent->toParent(_right);

                cout << "\n\ncomplex\n"
                     << "cur " << _cur
                     << "\nleft ** " << _left
                     << "\nright ** " << _right
                     << "\ntp ** " << _curParent
                     << "\nlp ** " << _leftParent
                     << "\nrp ** " << _rightParent << endl;
        */
        /*
            if (_cur->getLeft()->getTopSegment()->getArrayIndex() == 1219)
            {
              GappedBottomSegmentIteratorPtr bt = _leftParent->clone();
              cout << "bt copy " << bt << endl;
              bt->toRight();
              cout << "bt right " << bt << endl;
              _leftParent->toRight();
              cout << "\n after lp to right ** " << _leftParent;
              _rightParent->toLeft();
              cout << "\n after rp to left ** " << _rightParent << endl;
              scanDeletionCycle(topSegment);
        //      exit(10);
            }

              }
            }
        */
        _id = Complex;
    }

    return true;
}

bool Rearrangement::identifyDeletionFromLeftBreakpoint(TopSegmentIteratorPtr topSegment) {
    if (scanDeletionCycle(topSegment) == true && _leftParent->hasChild() == false) {
        _id = Deletion;
        return true;
    }
    _id = Invalid;
    return false;
}

pair<hal_index_t, hal_index_t> Rearrangement::getDeletedRange() const {
    pair<hal_index_t, hal_index_t> range;
    if (_leftParent->getReversed() == false) {
        range.first = _leftParent->getLeft()->getStartPosition();
        range.second = _leftParent->getRight()->getEndPosition();
    } else {
        range.first = _leftParent->getRight()->getEndPosition();
        range.second = _leftParent->getLeft()->getStartPosition();
    }
    assert(range.first <= range.second);
    assert(range.second - range.first == (hal_index_t)(getLength() - 1));
    return range;
}

bool Rearrangement::identifyInsertionFromLeftBreakpoint(TopSegmentIteratorPtr topSegment) {
    if (scanInsertionCycle(topSegment) == true && _cur->hasParent() == false) {
        _id = Insertion;
        return true;
    }
    _id = Invalid;
    return false;
}

pair<hal_index_t, hal_index_t> Rearrangement::getInsertedRange() const {
    pair<hal_index_t, hal_index_t> range;
    range.first = _cur->getLeft()->getStartPosition();
    range.second = _cur->getRight()->getStartPosition() + (hal_index_t)(_cur->getRight()->getLength() - 1);
    assert(_cur->getLeft()->tseg()->hasParent() == false);
    assert(_cur->getRight()->tseg()->hasParent() == false);
    if (range.first >= range.second) {
        swap(range.first, range.second);
    }
    assert(range.second - range.first == (hal_index_t)(getLength() - 1));
    return range;
}

pair<hal_index_t, hal_index_t> Rearrangement::getDuplicatedRange() const {
    pair<hal_index_t, hal_index_t> range;
    assert(_cur->hasParent() == true);
    _curParent->toParent(_cur);
    if (_curParent->getReversed() == false) {
        range.first = _curParent->getLeft()->getStartPosition();
        range.second = _curParent->getRight()->getEndPosition();
    } else {
        range.first = _curParent->getRight()->getEndPosition();
        range.second = _curParent->getLeft()->getStartPosition();
    }
    assert(range.first <= range.second);
    assert(range.second - range.first == (hal_index_t)(_curParent->getLength() - 1));
    return range;
}

bool Rearrangement::identifyNext() {
    assert(_cur->getReversed() == false);
    // don't like this.  need to refactor interface to make better
    // use of gapped iterators
    _top->copy(_cur->getRight());
    if (_top->getTopSegment()->getArrayIndex() < (hal_index_t)(_genome->getNumTopSegments() - 1)) {
        _top->toRight();
        identifyFromLeftBreakpoint(_top);
        return true;
    } else {
        return false;
    }
}

hal_size_t Rearrangement::getGapLengthThreshold() const {
    return _gapThreshold;
}

void Rearrangement::setGapLengthThreshold(hal_size_t threshold) {
    _gapThreshold = threshold;
    _cur = _genome->getGappedTopSegmentIterator(0, _gapThreshold, _atomic);
    _next = _cur->clone();
    _left = _cur->clone();
    _right = _left->clone();
    _leftParent = _parent->getGappedBottomSegmentIterator(0, 0, _gapThreshold, _atomic);
    _rightParent = _leftParent->clone();
    _curParent = _leftParent->clone();
    _nextParent = _leftParent->clone();
    _top = _cur->getLeft()->clone();
}

bool Rearrangement::getAtomic() const {
    return _atomic;
}

void Rearrangement::setAtomic(bool atomic) {
    _atomic = atomic;
    setGapLengthThreshold(0);
}

void Rearrangement::setNThreshold(double nThreshold) {
    if (nThreshold < 0 || nThreshold > 1) {
        throw hal_exception("Invalid nThreshold, " + std::to_string(nThreshold) + ", specified. "
                                                                                  "must be between 0 and 1.0");
    }
    _nThreshold = nThreshold;
}

double Rearrangement::getNThreshold() const {
    return _nThreshold;
}

void Rearrangement::resetStatus(TopSegmentIteratorPtr topSegment) {
    _id = Invalid;
    assert(topSegment.get());
    _genome = topSegment->getTopSegment()->getGenome();
    _parent = _genome->getParent();
    assert(_parent != NULL);

    _cur->setLeft(topSegment);
    _next->copy(_cur);
    _left->copy(_cur);
    _right->copy(_left);
    assert(_cur->getGapThreshold() == _gapThreshold);
    assert(_next->getGapThreshold() == _gapThreshold);
    assert(_left->getGapThreshold() == _gapThreshold);
    assert(_right->getGapThreshold() == _gapThreshold);
    assert(_leftParent->getGapThreshold() == _gapThreshold);
    assert(_rightParent->getGapThreshold() == _gapThreshold);
    assert(_curParent->getGapThreshold() == _gapThreshold);
    assert(_cur->getAtomic() == _atomic);
    assert(_next->getAtomic() == _atomic);
    assert(_left->getAtomic() == _atomic);
    assert(_right->getAtomic() == _atomic);
    assert(_leftParent->getAtomic() == _atomic);
    assert(_rightParent->getAtomic() == _atomic);
    assert(_curParent->getAtomic() == _atomic);
}

// Segment corresponds to no rearrangemnt.  This will happen when
// there is a rearrangement in the homolgous segment in its sibling
// genome.  In general, we can expect about half of segments to correspond
// to such cases.
bool Rearrangement::scanNothingCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);
    bool first = _cur->isFirst();
    bool last = _cur->isLast();

    if (_cur->hasParent() == false) {
        return false;
    }
    _curParent->toParent(_cur);
    if (first == false) {
        _left->toLeft();
        if (_left->hasParent() == false) {
            return false;
        }
        _leftParent->toParent(_left);
        if (_leftParent->adjacentTo(_curParent) == false) {
            return false;
        }
        if (_left->getParentReversed() == true) {
            if (_cur->getParentReversed() == false || _leftParent->rightOf(_curParent->getStartPosition()) == false) {
                return false;
            }
        } else {
            if (_cur->getParentReversed() == true || _leftParent->leftOf(_curParent->getStartPosition()) == false) {
                return false;
            }
        }
    }
    if (last == false) {
        _right->toRight();
        if (_right->hasParent() == false) {
            return false;
        }
        _rightParent->toParent(_right);
        if (_rightParent->adjacentTo(_curParent) == false) {
            return false;
        }
        if (_right->getParentReversed() == true) {
            if (_cur->getParentReversed() == false || _rightParent->leftOf(_curParent->getStartPosition()) == false) {
                return false;
            }
        } else {
            if (_cur->getParentReversed() == true || _rightParent->rightOf(_curParent->getStartPosition()) == false) {
                return false;
            }
        }
    }
    return last && first ? _cur->getParentReversed() : true;
}

// Segment is an reversed descendant of another Segment but
// otherwise no rearrangement.
bool Rearrangement::scanInversionCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);
    bool first = _cur->isFirst();
    bool last = _cur->isLast();

    if (_cur->hasParent() == false) {
        return false;
    }
    _curParent->toParent(_cur);
    if (first == false) {
        _left->toLeft();
        if (_left->hasParent() == false) {
            return false;
        }
        _leftParent->toParent(_left);
        if (_leftParent->adjacentTo(_curParent) == false) {
            return false;
        }
    }
    if (last == false) {
        _right->toRight();
        if (_right->hasParent() == false) {
            return false;
        }
        _rightParent->toParent(_right);
        if (_rightParent->adjacentTo(_curParent) == false) {
            return false;
        }
    }
    return _cur->getParentReversed();
}

// If true, _cur will store the insertion 'candidate'
// It must be further verified that this segment has no parent to
// distinguish between destination of transposition and insertion.
bool Rearrangement::scanInsertionCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);

    // eat up any adjacent insertions so they don't get double counted
    while (_next->hasParent() == false && _next->isLast() == false) {
        _right->copy(_next);
        _right->toRight();
        if (_right->hasParent() == false) {
            _next->copy(_right);
        } else {
            break;
        }
    }
    _right->copy(_next);
    assert(_next->equals(_cur) || _next->hasParent() == false);

    bool first = _cur->isFirst();
    bool last = _right->isLast();
    if (first && last) {
        return false;
    }

    // Case 1a) current segment is left endpoint.  we consider insertion
    // if right neighbour has parent
    if (first) {
        _right->toRight();
        if (_cur->hasParent() == false) {
            return true;
        } else if (_right->hasParent()) {
            _curParent->toParent(_cur);
            _rightParent->toParent(_right);
            return _rightParent->adjacentTo(_curParent) == false;
        }
    }

    // Case 1b) current segment is right endpoint.  we consider insertion
    // if left neighbour has parent
    else if (last) {
        _left->toLeft();
        if (_cur->hasParent() == false) {
            return true;
        } else if (_left->hasParent()) {
            _curParent->toParent(_cur);
            _leftParent->toParent(_left);
            return _leftParent->adjacentTo(_curParent) == false;
        }
    }

    // Case 2) current segment has a left neigbhour and a right neigbour
    else {
        _left->toLeft();
        _right->toRight();
        if (_left->hasParent() == true && _right->hasParent() == true) {
            _leftParent->toParent(_left);
            _rightParent->toParent(_right);
            // Case 2a) Parents are adjacent
            if (_leftParent->adjacentTo(_rightParent)) {
                return true;
            }
            // Case 2b) Left parent is endpoint
            else if (_leftParent->isFirst() || _leftParent->isLast()) {
                return _leftParent->getSequence() == _rightParent->getSequence();
            }

            // Case 2c) Right parent is endpoint
            else if (_rightParent->isFirst() || _rightParent->isLast()) {
                return _leftParent->getSequence() == _rightParent->getSequence();
            }
        }
    }

    return false;
}

// If true, _leftParent will store the deletion 'candidate'
// It must be further verified that this segment has no child to
// distinguish between source of transposition and deletion.
bool Rearrangement::scanDeletionCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);
    assert(_atomic != true || _cur->getNumSegments() == 1);

    bool first = _cur->isFirst();
    bool last = _cur->isLast();
    if (_cur->hasParent() == false || (first && last)) {
        return false;
    }

    // Case 1) current segment is a right endpoint.  we consider delection
    // if parent has neighbour
    if (last) {
        _leftParent->toParent(_cur);
        if (_leftParent->isFirst() == false) {
            _leftParent->toLeft();
            return true;
        }
        if (_leftParent->isLast() == false) {
            _leftParent->toRight();
            return true;
        }
    }

    // Case 2) Try to find deletion cycle by going right-up-left-left-down
    else {
        _leftParent->toParent(_cur);
        _right->toRight();

        assert(_leftParent->getGapThreshold() == _gapThreshold);
        assert(_cur->getGapThreshold() == _gapThreshold);
        assert(_atomic != true || _leftParent->getNumSegments() == 1);
        assert(_atomic != true || _right->getNumSegments() == 1);
        if (_right->hasParent() == false) {
            return false;
        }
        _rightParent->toParent(_right);

        if (_leftParent->getSequence() == _rightParent->getSequence()) {
            // don't care about inversions
            // so we make sure left is left of right and they are both positive
            if (_leftParent->getReversed() == true) {
                _leftParent->toReverse();
            }
            if (_rightParent->getReversed() == true) {
                _rightParent->toReverse();
            }
            if (_rightParent->getLeftArrayIndex() < _leftParent->getLeftArrayIndex()) {
                swap(_leftParent, _rightParent);
            }

            if (_leftParent->isLast()) {
                return false;
            }

            _leftParent->toRight();
            return _leftParent->adjacentTo(_rightParent);
        }
    }

    return false;
}

// If true, _leftParent will store the swapped segment (and _cur will store)
// the other half
// NEED TO REVISE WITH STRONGER CRITERIA -- right now any operation
// next to an endpoint can get confused with a translocation.
bool Rearrangement::scanTranslocationCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);
    bool first = _cur->isFirst();
    bool last = _cur->isLast();
    if (_cur->hasParent() == false || (!first && !last)) {
        return false;
    }

    _leftParent->toParent(_cur);
    bool pFirst = _leftParent->isFirst();
    // bool pLast = _leftParent->isLast();
    _rightParent->copy(_leftParent);

    first ? _right->toRight() : _right->toLeft();
    pFirst ? _rightParent->toRight() : _rightParent->toLeft();

    if (_right->hasParent() == false) {
        return true;
    } else {
        _curParent->toParent(_right);
        return _curParent->equals(_rightParent);
    }
    return false;
}

// leaves duplication on _cur and _right
bool Rearrangement::scanDuplicationCycle(TopSegmentIteratorPtr topSegment) {
    assert(topSegment.get());
    resetStatus(topSegment);
    return _cur->hasNextParalogy() == true && _cur->isCanonicalParalog() == false;
}
