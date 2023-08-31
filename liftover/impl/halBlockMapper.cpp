/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halBlockMapper.h"
#include "hal.h"
#include "halSegmentMapper.h"
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace hal;
using namespace std;

hal_size_t BlockMapper::_maxAdjScan = 1;

BlockMapper::BlockMapper() {
}

BlockMapper::~BlockMapper() {
    erase();
}

void BlockMapper::erase() {
    _segSet.clear();
    _adjSet.clear();
    _downwardPath.clear();
    _upwardPath.clear();
}

void BlockMapper::init(const Genome *refGenome, const Genome *queryGenome, hal_index_t absRefFirst, hal_index_t absRefLast,
                       bool targetReversed, bool doDupes, hal_size_t minLength, bool mapTargetAdjacencies,
                       const Genome *coalescenceLimit) {
    erase();
    _absRefFirst = absRefFirst;
    _absRefLast = absRefLast;
    _targetReversed = targetReversed;
    _doDupes = doDupes;
    _minLength = minLength;
    _mapAdj = mapTargetAdjacencies;
    _refGenome = refGenome;
    _refSequence = refGenome->getSequenceBySite(_absRefFirst);
    assert(_refSequence == refGenome->getSequenceBySite(_absRefLast));
    _queryGenome = queryGenome;

    set<const Genome *> inputSet;
    inputSet.insert(_refGenome);
    inputSet.insert(_queryGenome);
    _mrca = getLowestCommonAncestor(inputSet);

    if (coalescenceLimit == NULL) {
        _coalescenceLimit = _mrca;
    } else {
        _coalescenceLimit = coalescenceLimit;
    }

    // The path between the coalescence limit (the highest point in the
    // tree) and the query genome is needed to traverse down into the
    // correct children.
    inputSet.clear();
    inputSet.insert(_queryGenome);
    inputSet.insert(_coalescenceLimit);
    getGenomesInSpanningTree(inputSet, _downwardPath);

    // similarly, the upward path is needed to get the adjacencies properly.
    inputSet.clear();
    inputSet.insert(_refGenome);
    inputSet.insert(_coalescenceLimit);
    getGenomesInSpanningTree(inputSet, _upwardPath);
}

void BlockMapper::map() {
    SegmentIteratorPtr refSeg;
    hal_index_t lastIndex;
    if ((_mrca == _refGenome) && (_refGenome != _queryGenome)) {
        refSeg = _refGenome->getBottomSegmentIterator();
        lastIndex = _refGenome->getNumBottomSegments();
    } else {
        refSeg = _refGenome->getTopSegmentIterator();
        lastIndex = _refGenome->getNumTopSegments();
    }

    refSeg->toSite(_absRefFirst, false);
    hal_offset_t startOffset = _absRefFirst - refSeg->getStartPosition();
    hal_offset_t endOffset = 0;
    if (_absRefLast <= refSeg->getEndPosition()) {
        endOffset = refSeg->getEndPosition() - _absRefLast;
    }
    refSeg->slice(startOffset, endOffset);

    assert(refSeg->getStartPosition() == _absRefFirst);
    assert(refSeg->getEndPosition() <= _absRefLast);

    while (refSeg->getArrayIndex() < lastIndex && refSeg->getStartPosition() <= _absRefLast) {
        if (_targetReversed == true) {
            refSeg->toReverseInPlace();
        }
        halMapSegment(refSeg.get(), _segSet, _queryGenome, &_downwardPath, _doDupes, _minLength, _coalescenceLimit, _mrca);
        if (_targetReversed == true) {
            refSeg->toReverseInPlace();
        }
        refSeg->toRight(_absRefLast);
    }

    if (_mapAdj) {
        assert(_targetReversed == false);
        MappedSegmentSet::const_iterator i;
        for (i = _segSet.begin(); i != _segSet.end(); ++i) {
            if (_adjSet.find(*i) == _adjSet.end()) {
                mapAdjacencies(i);
            }
        }
    }
}

void BlockMapper::mapAdjacencies(MappedSegmentSet::const_iterator segIt) {
    assert(_segSet.empty() == false && segIt != _segSet.end());
    MappedSegmentPtr mappedQuerySeg(*segIt);
    hal_index_t maxIndex;
    hal_index_t minIndex;
    SegmentIteratorPtr queryIt = makeIterator(mappedQuerySeg, minIndex, maxIndex);
    MappedSegmentSet backResults;
    MappedSegmentSet::const_iterator segNext = segIt;
    if (queryIt->getReversed()) {
        segNext = segNext == _segSet.begin() ? _segSet.end() : --segNext;
    } else {
        ++segNext;
    }

    hal_size_t iter = 0;
    queryIt->toRight();
    while (queryIt->getArrayIndex() >= minIndex && queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan) {
        bool wasCut = false;
        if (segNext != _segSet.end()) {
            // if cut returns nothing, then the region in question is covered
            // by segNext (ie already mapped).
            wasCut = cutByNext(queryIt.get(), segNext->get()->getTarget(), !queryIt->getReversed());
        }
        if (wasCut == true) {
            break;
        }
        size_t backSize = backResults.size();
        assert(queryIt->getArrayIndex() >= 0);
        halMapSegment(queryIt.get(), backResults, _refGenome, &_upwardPath, _doDupes, _minLength);
        // something was found, that's good enough.
        if (backResults.size() > backSize) {
            break;
        }
        queryIt->toRight();
        ++iter;
    }

    queryIt = makeIterator(mappedQuerySeg, minIndex, maxIndex);

    MappedSegmentSet::const_iterator segPrev = segIt;
    if (queryIt->getReversed()) {
        ++segPrev;
    } else {
        segPrev = segPrev == _segSet.begin() ? _segSet.end() : --segPrev;
    }
    iter = 0;
    queryIt->toLeft();
    while (queryIt->getArrayIndex() >= minIndex && queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan) {
        bool wasCut = false;
        if (segPrev != _segSet.end()) {
            // if cut returns nothing, then the region in question is covered
            // by segPrev (ie already mapped).
            wasCut = cutByNext(queryIt.get(), segPrev->get()->getTarget(), queryIt->getReversed());
        }
        if (wasCut == true) {
            break;
        }
        size_t backSize = backResults.size();
        halMapSegment(queryIt.get(), backResults, _refGenome, &_upwardPath, _doDupes, _minLength);
        // something was found, that's good enough.
        if (backResults.size() > backSize) {
            break;
        }
        queryIt->toLeft();
        ++iter;
    }

    MappedSegmentSet outSet;
    // flip the results and copy back to our main set.
    for (MappedSegmentSet::iterator i = backResults.begin(); i != backResults.end(); ++i) {
        MappedSegmentPtr mseg(*i);
        if (mseg->getSequence() == _refSequence) {
            mseg->flip();
            const SlicedSegment *refSeg = mseg->getSource();
            if (refSeg->getReversed()) {
                mseg->fullReverse();
            }

            MappedSegmentSet::const_iterator j = _segSet.lower_bound(*i);
            bool overlaps = false;
            if (j != _segSet.begin()) {
                --j;
            }
            for (size_t count = 0; count < 3 && j != _segSet.end() && !overlaps; ++count, ++j) {
                overlaps = mseg->overlaps((*j)->getStartPosition()) || mseg->overlaps((*j)->getEndPosition()) ||
                           (*j)->overlaps(mseg->getStartPosition()) || (*j)->overlaps(mseg->getEndPosition());
            }
            if (overlaps == false) {
                outSet.insert(mseg);
            }
        }
    }

    // clean up dupes before adding to output
    for (MappedSegmentSet::iterator i = outSet.begin(); i != outSet.end();) {
        
        // find the equivalence class of identical target intervals
        MappedSegmentSet::iterator j = i;
        ++j;
        hal_index_t copies = 1;
        while (j != outSet.end() && ((*j)->getStartPosition() == (*i)->getStartPosition() ||
                                     (*j)->getEndPosition() == (*i)->getStartPosition())) {
            ++j;
            ++copies;
        }

        // choose the best copy based on Source distance to input (ie nearest in screen coordinates)
        MappedSegmentSet::iterator best;
        hal_index_t best_delta = numeric_limits<hal_index_t>::max();
        for (MappedSegmentSet::iterator k = i; k !=j; ++k) {
            hal_index_t delta = min(abs((*k)->getSource()->getStartPosition() - (*segIt)->getSource()->getStartPosition()),
                                    abs((*k)->getSource()->getEndPosition() - (*segIt)->getSource()->getStartPosition()));
            if (delta < best_delta) {
                best_delta = delta;
                best = k;
            }
        }

        _segSet.insert(*best);
        _adjSet.insert(*best);
       
        i = j;
    }

}

SegmentIteratorPtr BlockMapper::makeIterator(MappedSegmentPtr &mappedSegment, hal_index_t &minIndex, hal_index_t &maxIndex) {
    SegmentIteratorPtr segIt;
    if (mappedSegment->isTop()) {
        segIt = mappedSegment->getGenome()->getTopSegmentIterator(mappedSegment->getArrayIndex());
        minIndex = segIt->getSequence()->getTopSegmentArrayIndex();
        maxIndex = minIndex + (hal_index_t)segIt->getSequence()->getNumTopSegments();
    } else {
        segIt = mappedSegment->getGenome()->getBottomSegmentIterator(mappedSegment->getArrayIndex());
        minIndex = segIt->getSequence()->getBottomSegmentArrayIndex();
        maxIndex = minIndex + (hal_index_t)segIt->getSequence()->getNumBottomSegments();
    }

    if (mappedSegment->getReversed()) {
        segIt->toReverse();
    }
    segIt->slice(mappedSegment->getStartOffset(), mappedSegment->getEndOffset());

    assert(segIt->getGenome() == mappedSegment->getGenome());
    assert(segIt->getArrayIndex() == mappedSegment->getArrayIndex());
    assert(segIt->getStartOffset() == mappedSegment->getStartOffset());
    assert(segIt->getEndOffset() == mappedSegment->getEndOffset());
    assert(segIt->getReversed() == mappedSegment->getReversed());

    return segIt;
}

bool BlockMapper::cutByNext(SlicedSegment *query, const SlicedSegment *nextSeg, bool right) {
    assert(query->getGenome() == nextSeg->getGenome());
    assert(query->isTop() == nextSeg->isTop());
    bool wasCut = false;

    if (query->getArrayIndex() == nextSeg->getArrayIndex()) {
        hal_offset_t so1 = query->getStartOffset();
        hal_offset_t eo1 = query->getEndOffset();
        if (query->getReversed()) {
            swap(so1, eo1);
        }
        hal_offset_t so2 = nextSeg->getReversed() ? nextSeg->getEndOffset() : nextSeg->getStartOffset();

        if (right) {
            // overlap on start position.  we zap
            if (so1 >= so2) {
                wasCut = true;
            } else {
                hal_index_t e1 = max(query->getEndPosition(), query->getStartPosition());
                hal_index_t s2 = min(nextSeg->getEndPosition(), nextSeg->getStartPosition());
                // end position of query overlaps next seg.  so we cut it.
                if (e1 >= s2) {
                    hal_index_t delta = 1 + e1 - s2;
                    assert(delta < (hal_index_t)query->getLength());
                    hal_offset_t newEndOffset = eo1 + delta;
                    hal_offset_t newStartOffset = so1;
                    if (query->getReversed() == true) {
                        swap(newEndOffset, newStartOffset);
                    }
                    query->slice(newStartOffset, newEndOffset);
                }
            }
        } else {
            hal_index_t s1 = min(query->getEndPosition(), query->getStartPosition());
            hal_index_t e1 = max(query->getEndPosition(), query->getStartPosition());
            hal_index_t e2 = max(nextSeg->getEndPosition(), nextSeg->getStartPosition());

            // overlap on end position.  we zap
            if (e1 <= e2) {
                wasCut = true;
            } else {
                // end position of query overlaps next seg.  so we cut it.
                if (s1 <= e2) {
                    hal_index_t delta = 1 + e2 - s1;
                    assert(delta < (hal_index_t)query->getLength());
                    hal_offset_t newStartOffset = so1 + delta;
                    hal_offset_t newEndOffset = eo1;
                    if (query->getReversed() == true) {
                        swap(newEndOffset, newStartOffset);
                    }
                    query->slice(newStartOffset, newEndOffset);
                }
            }
        }
    }
    return wasCut;
}

void BlockMapper::extractSegment(MappedSegmentSet::iterator start, const MappedSegmentSet &paraSet,
                                 vector<MappedSegmentPtr> &fragments, MappedSegmentSet *startSet,
                                 const set<hal_index_t> &targetCutPoints, set<hal_index_t> &queryCutPoints) {
    fragments.clear();
    fragments.push_back(*start);
    const Sequence *startSeq = (*start)->getSequence();

    vector<MappedSegmentSet::iterator> vector1;
    vector<MappedSegmentSet::iterator> *v1 = &vector1;
    vector<MappedSegmentSet::iterator> vector2;
    vector<MappedSegmentSet::iterator> *v2 = &vector2;
    vector<MappedSegmentSet::iterator> toErase;

    v1->push_back(start);
    MappedSegmentSet::iterator next = start;
    ++next;

    // equivalence class based on start set
    while (next != startSet->end() && equalTargetStart(*v1->back(), *next)) {
        assert((*next)->getLength() == (*v1->back())->getLength());
        v1->push_back(next);
        ++next;
    }

    while (next != startSet->end()) {
        // equivalence class based on next element
        while (next != startSet->end() && (v2->empty() || equalTargetStart(*v2->back(), *next)) && v2->size() < v1->size()) {
            assert(v2->empty() || (*next)->getLength() == (*v2->back())->getLength());
            v2->push_back(next);
            ++next;
        }

        // check if all elements of each class are compatible for a merge
        assert(v1->size() > 0 && v2->size() > 0);
        bool canMerge = v1->size() == v2->size();
        for (size_t i = 0; i < v1->size() && canMerge; ++i) {
            canMerge = (v1->size() == v2->size() && (*v2->at(i))->getSequence() == startSeq &&
                        (*v1->at(i))->canMergeRightWith(*v2->at(i), &queryCutPoints, &targetCutPoints) &&
                        (paraSet.find(*v1->at(i)) == paraSet.end()) == (paraSet.find(*v2->at(i)) == paraSet.end()));
        }
        if (canMerge == true) {
            // add the next element and flag for deletion from start set
            fragments.push_back(*v2->at(0));
            toErase.push_back(v2->at(0));
        } else {
            break;
        }
        v1->clear();
        swap(v1, v2);
    }

    // we extracted a paralgous region along query coordinates.  we
    // add its right query coordinate to the cutset to make sure none
    // of the other paralogs ever get merged beyond it.
    if (v1->size() > 1) {
        queryCutPoints.insert(max(fragments.back()->getStartPosition(), fragments.back()->getEndPosition()));
    }

    for (size_t i = 0; i < toErase.size(); ++i) {
        startSet->erase(toErase[i]);
    }

    assert(fragments.front()->getSequence() == fragments.back()->getSequence());
}
