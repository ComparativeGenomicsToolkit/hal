/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halSegmentMapper.h"
#include "halBottomSegmentIterator.h"
#include "halCommon.h"
#include "halMappedSegment.h"
#include "halSegment.h"
#include "halSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace hal;

enum OverlapCat { Same, Disjoint, AContainsB, BContainsA, AOverlapsLeftOfB, BOverlapsLeftOfA };

static hal_size_t mapSelf(MappedSegmentPtr mappedSeg, list<MappedSegmentPtr> &results, hal_size_t minLength);

// note: takes smart pointer as it maybe added to the results
static hal_size_t mapUp(MappedSegmentPtr mappedSeg, list<MappedSegmentPtr> &results, bool doDupes, hal_size_t minLength) {
    const Genome *parent = mappedSeg->getGenome()->getParent();
    assert(parent != NULL);
    hal_size_t added = 0;
    if (mappedSeg->isTop() == true) {
        BottomSegmentIteratorPtr botSegIt = parent->getBottomSegmentIterator();
        TopSegmentIteratorPtr topSegIt = mappedSeg->targetAsTop();
        if (topSegIt->tseg()->hasParent() == true && topSegIt->getLength() >= minLength &&
            (doDupes == true || topSegIt->tseg()->isCanonicalParalog() == true)) {
            botSegIt->toParent(topSegIt);
            mappedSeg->setTarget(std::dynamic_pointer_cast<SegmentIterator>(botSegIt));
            results.push_back(mappedSeg);
            ++added;
        }
    } else {
        hal_index_t rightCutoff = mappedSeg->getEndPosition();
        BottomSegmentIteratorPtr botSegIt = mappedSeg->targetAsBottom();
        hal_index_t startOffset = (hal_index_t)botSegIt->getStartOffset();
        hal_index_t endOffset = (hal_index_t)botSegIt->getEndOffset();
        TopSegmentIteratorPtr topSegIt = mappedSeg->getGenome()->getTopSegmentIterator();
        topSegIt->toParseUp(botSegIt);
        do {
            TopSegmentIteratorPtr newTopSegIt = topSegIt->clone();

            // we map the new target back to see how the offsets have
            // changed.  these changes are then applied to the source segment
            // as deltas
            BottomSegmentIteratorPtr backBotSegIt = botSegIt->clone();
            backBotSegIt->toParseDown(newTopSegIt);
            hal_index_t startBack = (hal_index_t)backBotSegIt->getStartOffset();
            hal_index_t endBack = (hal_index_t)backBotSegIt->getEndOffset();
            assert(startBack >= startOffset);
            assert(endBack >= endOffset);
            SegmentIteratorPtr newSourceSegIt = mappedSeg->sourceClone();
            hal_index_t startDelta = startBack - startOffset;
            hal_index_t endDelta = endBack - endOffset;
            assert((hal_index_t)newSourceSegIt->getLength() > startDelta + endDelta);
            newSourceSegIt->slice(newSourceSegIt->getStartOffset() + startDelta, newSourceSegIt->getEndOffset() + endDelta);

            MappedSegmentPtr newMappedSeg(new MappedSegment(newSourceSegIt, newTopSegIt));

            assert(newMappedSeg->isTop() == true);
            assert(newMappedSeg->getSource()->getGenome() == mappedSeg->getSource()->getGenome());

            added += mapUp(newMappedSeg, results, doDupes, minLength);
            // stupid that we have to make this check but odn't want to
            // make fundamental api change now
            if (topSegIt->getEndPosition() != rightCutoff) {
                topSegIt->toRight(rightCutoff);
            } else {
                break;
            }
        } while (true);
    }
    return added;
}

// Map the input segments up until reaching the target genome. If the
// target genome is below the source genome, fail miserably.
// Destructive to any data in the input or results list.
static hal_size_t mapRecursiveUp(list<MappedSegmentPtr> &input, list<MappedSegmentPtr> &results, const Genome *tgtGenome,
                                 hal_size_t minLength) {
    list<MappedSegmentPtr> *inputPtr = &input;
    list<MappedSegmentPtr> *outputPtr = &results;

    if (inputPtr->empty() || (*inputPtr->begin())->getGenome() == tgtGenome) {
        results = *inputPtr;
        return 0;
    }

    const Genome *curGenome = (*inputPtr->begin())->getGenome();
    assert(curGenome != NULL);
    const Genome *nextGenome = curGenome->getParent();

    if (nextGenome == NULL) {
        throw hal_exception("Reached top of tree when attempting to recursively map up from " + curGenome->getName() + " to " +
                            tgtGenome->getName());
    }

    // Map all segments to the parent.
    list<MappedSegmentPtr>::iterator i = inputPtr->begin();
    for (; i != inputPtr->end(); ++i) {
        assert((*i)->getGenome() == curGenome);
        mapUp(*i, *outputPtr, true, minLength);
    }

    if (nextGenome != tgtGenome) {
        // Continue the recursion.
        swap(inputPtr, outputPtr);
        outputPtr->clear();
        mapRecursiveUp(*inputPtr, *outputPtr, tgtGenome, minLength);
    }

    if (outputPtr != &results) {
        results = *outputPtr;
    }

    results.sort(MappedSegment::LessSourcePtr());
    results.unique(MappedSegment::EqualToPtr());
    return results.size();
}

// note: takes smart pointer as it maybe added to the results
static hal_size_t mapDown(MappedSegmentPtr mappedSeg, hal_size_t childIndex, list<MappedSegmentPtr> &results,
                          hal_size_t minLength) {
    const Genome *child = mappedSeg->getGenome()->getChild(childIndex);
    assert(child != NULL);
    hal_size_t added = 0;
    if (mappedSeg->isTop() == false) {
        TopSegmentIteratorPtr topSegIt = child->getTopSegmentIterator();
        SegmentIteratorPtr targetSegIt = mappedSeg->getTargetIteratorPtr();
        BottomSegmentIteratorPtr botSegIt = std::dynamic_pointer_cast<BottomSegmentIterator>(targetSegIt);

        if (botSegIt->bseg()->hasChild(childIndex) == true && botSegIt->getLength() >= minLength) {
            topSegIt->toChild(botSegIt, childIndex);
            mappedSeg->setTarget(std::dynamic_pointer_cast<SegmentIterator>(topSegIt));
            results.push_back(MappedSegmentPtr(mappedSeg));
            ++added;
        }
    } else {
        hal_index_t rightCutoff = mappedSeg->getEndPosition();
        TopSegmentIteratorPtr topSegIt = mappedSeg->targetAsTop();
        hal_index_t startOffset = (hal_index_t)topSegIt->getStartOffset();
        hal_index_t endOffset = (hal_index_t)topSegIt->getEndOffset();
        BottomSegmentIteratorPtr botSegIt = mappedSeg->getGenome()->getBottomSegmentIterator();
        botSegIt->toParseDown(topSegIt);
        do {
            BottomSegmentIteratorPtr newBotSegIt = botSegIt->clone();

            // we map the new target back to see how the offsets have
            // changed.  these changes are then applied to the source segment
            // as deltas
            TopSegmentIteratorPtr backTopSegIt = topSegIt->clone();
            backTopSegIt->toParseUp(newBotSegIt);
            hal_index_t startBack = (hal_index_t)backTopSegIt->getStartOffset();
            hal_index_t endBack = (hal_index_t)backTopSegIt->getEndOffset();
            assert(startBack >= startOffset);
            assert(endBack >= endOffset);
            SegmentIteratorPtr newSourceSegIt = mappedSeg->sourceClone();
            hal_index_t startDelta = startBack - startOffset;
            hal_index_t endDelta = endBack - endOffset;
            assert((hal_index_t)newSourceSegIt->getLength() > startDelta + endDelta);
            newSourceSegIt->slice(newSourceSegIt->getStartOffset() + startDelta, newSourceSegIt->getEndOffset() + endDelta);

            MappedSegmentPtr newMappedSeg(new MappedSegment(newSourceSegIt, newBotSegIt));

            assert(newMappedSeg->isTop() == false);
            assert(newMappedSeg->getSource()->getGenome() == mappedSeg->getSource()->getGenome());

            added += mapDown(newMappedSeg, childIndex, results, minLength);

            // stupid that we have to make this check but odn't want to
            // make fundamental api change now
            if (botSegIt->getEndPosition() != rightCutoff) {
                botSegIt->toRight(rightCutoff);
            } else {
                break;
            }
        } while (true);
    }
    return added;
}

// Map the input segments down until reaching the target genome. If the
// target genome is above the source genome, fail miserably.
// Destructive to any data in the input or results list.
static hal_size_t mapRecursiveDown(list<MappedSegmentPtr> &input, list<MappedSegmentPtr> &results, const Genome *tgtGenome,
                                   const set<string> &namesOnPath, bool doDupes, hal_size_t minLength) {
    list<MappedSegmentPtr> *inputPtr = &input;
    list<MappedSegmentPtr> *outputPtr = &results;

    if (inputPtr->empty()) {
        results = *inputPtr;
        return 0;
    }

    const Genome *curGenome = (*inputPtr->begin())->getGenome();
    assert(curGenome != NULL);
    if (curGenome == tgtGenome) {
        results = *inputPtr;
        return 0;
    }

    // Find the correct child to move down into.
    const Genome *nextGenome = NULL;
    hal_size_t nextChildIndex = numeric_limits<hal_size_t>::max();
    const Alignment *alignment = curGenome->getAlignment();
    vector<string> childNames = alignment->getChildNames(curGenome->getName());
    for (hal_size_t child = 0; nextGenome == NULL && child < childNames.size(); ++child) {
        if (childNames[child] == tgtGenome->getName() || namesOnPath.find(childNames[child]) != namesOnPath.end()) {
            const Genome *childGenome = curGenome->getChild(child);
            nextGenome = childGenome;
            nextChildIndex = child;
        }
    }

    if (nextGenome == NULL) {
        throw hal_exception("Could not find correct child that leads from " + curGenome->getName() + " to " +
                            tgtGenome->getName());
    }

    assert(nextGenome->getParent() == curGenome);

    // Map the actual segments down.
    list<MappedSegmentPtr>::iterator i = inputPtr->begin();
    for (; i != inputPtr->end(); ++i) {
        assert((*i)->getGenome() == curGenome);
        mapDown(*i, nextChildIndex, *outputPtr, minLength);
    }

    // Find paralogs.
    if (doDupes == true) {
        swap(inputPtr, outputPtr);
        outputPtr->clear();
        list<MappedSegmentPtr>::iterator i = inputPtr->begin();
        for (; i != inputPtr->end(); ++i) {
            assert((*i)->getGenome() == nextGenome);
            mapSelf(*i, *outputPtr, minLength);
        }
    }

    if (nextGenome != tgtGenome) {
        // Continue the recursion.
        swap(inputPtr, outputPtr);
        outputPtr->clear();
        mapRecursiveDown(*inputPtr, *outputPtr, tgtGenome, namesOnPath, doDupes, minLength);
    }

    if (outputPtr != &results) {
        results = *outputPtr;
    }

    results.sort(MappedSegment::LessSourcePtr());
    results.unique(MappedSegment::EqualToPtr());
    return results.size();
}

// note: takes smart pointer as it maybe added to the results
static hal_size_t mapSelf(MappedSegmentPtr mappedSeg, list<MappedSegmentPtr> &results, hal_size_t minLength) {
    hal_size_t added = 0;
    if (mappedSeg->isTop() == true) {
        SegmentIteratorPtr target = mappedSeg->getTargetIteratorPtr();
        SegmentIteratorPtr source = mappedSeg->getSourceIteratorPtr();
        TopSegmentIteratorPtr top = std::dynamic_pointer_cast<TopSegmentIterator>(target);
        TopSegmentIteratorPtr topCopy = top->clone();
        do {
            // FIXME: why isn't clone() polymorphic?
            SegmentIteratorPtr newSource;
            if (source->isTop()) {
                newSource = std::dynamic_pointer_cast<TopSegmentIterator>(source)->clone();
            } else {
                newSource = std::dynamic_pointer_cast<BottomSegmentIterator>(source)->clone();
            }
            TopSegmentIteratorPtr newTop = topCopy->clone();
            MappedSegmentPtr newMappedSeg(new MappedSegment(newSource, newTop));
            assert(newMappedSeg->getGenome() == mappedSeg->getGenome());
            assert(newMappedSeg->getSource()->getGenome() == mappedSeg->getSource()->getGenome());
            results.push_back(newMappedSeg);
            ++added;
            if (topCopy->tseg()->hasNextParalogy()) {
                topCopy->toNextParalogy();
            }
        } while (topCopy->tseg()->hasNextParalogy() == true && topCopy->getLength() >= minLength &&
                 topCopy->getArrayIndex() != top->getArrayIndex());
    } else if (mappedSeg->getGenome()->getParent() != NULL) {
        hal_index_t rightCutoff = mappedSeg->getEndPosition();
        BottomSegmentIteratorPtr bottom = mappedSeg->targetAsBottom();
        hal_index_t startOffset = (hal_index_t)bottom->getStartOffset();
        hal_index_t endOffset = (hal_index_t)bottom->getEndOffset();
        TopSegmentIteratorPtr top = mappedSeg->getGenome()->getTopSegmentIterator();
        top->toParseUp(bottom);
        do {
            TopSegmentIteratorPtr topNew = top->clone();

            // we map the new target back to see how the offsets have
            // changed.  these changes are then applied to the source segment
            // as deltas
            BottomSegmentIteratorPtr bottomBack = bottom->clone();
            bottomBack->toParseDown(topNew);
            hal_index_t startBack = (hal_index_t)bottomBack->getStartOffset();
            hal_index_t endBack = (hal_index_t)bottomBack->getEndOffset();
            assert(startBack >= startOffset);
            assert(endBack >= endOffset);
            SegmentIteratorPtr newSource = mappedSeg->sourceClone();
            hal_index_t startDelta = startBack - startOffset;
            hal_index_t endDelta = endBack - endOffset;
            assert((hal_index_t)newSource->getLength() > startDelta + endDelta);
            newSource->slice(newSource->getStartOffset() + startDelta, newSource->getEndOffset() + endDelta);

            MappedSegmentPtr newMappedSeg(new MappedSegment(newSource, topNew));

            assert(newMappedSeg->isTop() == true);
            assert(newMappedSeg->getSource()->getGenome() == mappedSeg->getSource()->getGenome());

            added += mapSelf(newMappedSeg, results, minLength);
            // stupid that we have to make this check but odn't want to
            // make fundamental api change now
            if (top->getEndPosition() != rightCutoff) {
                top->toRight(rightCutoff);
            } else {
                break;
            }
        } while (true);
    }
    return added;
}

static OverlapCat slowOverlap(const SlicedSegment *sA, const SlicedSegment *sB) {
    hal_index_t startA = sA->getStartPosition();
    hal_index_t endA = sA->getEndPosition();
    hal_index_t startB = sB->getStartPosition();
    hal_index_t endB = sB->getEndPosition();
    if (startA > endA) {
        swap(startA, endA);
    }
    if (startB > endB) {
        swap(startB, endB);
    }

    if (endA < startB || startA > endB) {
        return Disjoint;
    } else if (startA == startB && endA == endB) {
        return Same;
    } else if (startA >= startB && endA <= endB) {
        return BContainsA;
    } else if (startB >= startA && endB <= endA) {
        return AContainsB;
    } else if (startA <= startB && endA < endB) {
        return AOverlapsLeftOfB;
    }
    assert(startB <= startA && endB < endA);
    return BOverlapsLeftOfA;
}

static void getOverlapBounds(MappedSegmentPtr &seg, MappedSegmentSet &results, MappedSegmentSet::iterator &leftBound,
                             MappedSegmentSet::iterator &rightBound) {
    if (results.size() <= 2) {
        leftBound = results.begin();
        rightBound = results.end();
    } else {
        MappedSegmentSet::iterator i = results.lower_bound(seg);
        leftBound = i;
        if (leftBound != results.begin()) {
            --leftBound;
        }
        MappedSegmentSet::iterator iprev;
        MappedSegmentSet::key_compare resLess = results.key_comp();
        while (leftBound != results.begin()) {
            iprev = leftBound;
            --iprev;
            if (leftBound == results.end() || !resLess(*iprev, *leftBound)) {
                leftBound = iprev;
            } else {
                break;
            }
        }
        for (; leftBound != results.begin(); --leftBound) {
            if (leftBound != results.end() && slowOverlap(seg.get()->getTarget(), leftBound->get()->getTarget()) == Disjoint) {
                break;
            }
        }
        rightBound = i;
        if (rightBound != results.end()) {
            for (++rightBound; rightBound != results.end(); ++rightBound) {
                if (slowOverlap(seg.get()->getTarget(), rightBound->get()->getTarget()) == Disjoint) {
                    break;
                }
            }
        }
    }
}

static void clipAagainstB(MappedSegmentPtr segA, MappedSegmentPtr segB, OverlapCat overlapCat,
                          vector<MappedSegmentPtr> &clippedSegs) {
    assert(overlapCat != Same && overlapCat != Disjoint && overlapCat != BContainsA);

    hal_index_t startA = segA->getStartPosition();
    hal_index_t endA = segA->getEndPosition();
    hal_index_t startB = segB->getStartPosition();
    hal_index_t endB = segB->getEndPosition();
    if (startA > endA) {
        swap(startA, endA);
    }
    if (startB > endB) {
        swap(startB, endB);
    }
    MappedSegmentPtr left = segA;
    MappedSegmentPtr middle = MappedSegmentPtr(segA->clone());
    MappedSegmentPtr right;

    hal_index_t startO = segA->getStartOffset();
    hal_index_t endO = segA->getEndOffset();
    hal_index_t length = segA->getLength();
    hal_index_t leftSize = std::max((hal_index_t)0, startB - startA);
    hal_index_t rightSize = std::max((hal_index_t)0, endA - endB);
    hal_index_t middleSize = length - leftSize - rightSize;
    if (rightSize > 0) {
        right = MappedSegmentPtr(segA->clone());
    }

    assert(overlapCat == AOverlapsLeftOfB || overlapCat == BOverlapsLeftOfA || overlapCat == AContainsB);

    assert(leftSize >= 0 && rightSize >= 0 && middleSize >= 0);
    hal_index_t leftSlice = 0;
    hal_index_t rightSlice = 0;
    if (leftSize > 0) {
        leftSlice = 0;
        rightSlice = (length - leftSize);
        if (left->getReversed()) {
            swap(leftSlice, rightSlice);
        }
        left->slice(startO + leftSlice, endO + rightSlice);
        assert(left->getLength() == (hal_size_t)leftSize);
        assert(min(left->getStartPosition(), left->getEndPosition()) == startA);
        assert(max(left->getStartPosition(), left->getEndPosition()) == startB - 1);
    } else {
        middle = segA;
    }

    leftSlice = leftSize;
    rightSlice = rightSize;
    if (middle->getReversed()) {
        swap(leftSlice, rightSlice);
    }
    middle->slice(startO + leftSlice, endO + rightSlice);
    assert(middle->getLength() == (hal_size_t)middleSize);
    assert(min(middle->getStartPosition(), middle->getEndPosition()) == max(startB, startA));
    assert(max(middle->getStartPosition(), middle->getEndPosition()) == min(endB, endA));
    if (middle.get() != segA.get()) {
        assert(leftSize > 0);
        assert(middle->getLength() == middle->getSource()->getLength());
        clippedSegs.push_back(middle);
    }

    if (rightSize > 0) {
        leftSlice = leftSize + middleSize;
        rightSlice = 0;
        if (right->getReversed()) {
            swap(leftSlice, rightSlice);
        }
        right->slice(startO + leftSlice, endO + rightSlice);
        assert(right->getLength() == (hal_size_t)rightSize);
        assert(min(right->getStartPosition(), right->getEndPosition()) == endB + 1);
        assert(max(right->getStartPosition(), right->getEndPosition()) == endA);
        assert(right->getLength() == right->getSource()->getLength());
        clippedSegs.push_back(right);
    }
    assert(segA->getLength() == segA->getSource()->getLength());
}

static void insertAndBreakOverlaps(MappedSegmentPtr seg, MappedSegmentSet &results) {
    assert(seg->getLength() == seg->getSource()->getLength());
    list<MappedSegmentPtr> inputSegs;
    vector<MappedSegmentPtr> clippedSegs;

    // 1) compute invariant range in set of candidate overalaps
    MappedSegmentSet::iterator leftBound;
    MappedSegmentSet::iterator rightBound;
    getOverlapBounds(seg, results, leftBound, rightBound);
    bool leftBegin = leftBound == results.begin();

    // 2) cut seg by each segment in range
    list<MappedSegmentPtr>::iterator inputIt;
    OverlapCat oc;
    inputSegs.push_back(seg);
    MappedSegmentSet::iterator resIt;
    for (resIt = leftBound; resIt != rightBound; ++resIt) {
        for (inputIt = inputSegs.begin(); inputIt != inputSegs.end(); ++inputIt) {
            oc = slowOverlap(inputIt->get()->getTarget(), resIt->get()->getTarget());
            if (oc == AContainsB || oc == AOverlapsLeftOfB || oc == BOverlapsLeftOfA) {
                clippedSegs.clear();
                clipAagainstB(*inputIt, *resIt, oc, clippedSegs);
                inputSegs.insert(inputSegs.end(), clippedSegs.begin(), clippedSegs.end());
            }
        }
    }

    // 3) cut results by input list
    for (inputIt = inputSegs.begin(); inputIt != inputSegs.end(); ++inputIt) {
        assert((*inputIt)->getLength() == (*inputIt)->getSource()->getLength());
        resIt = leftBegin ? results.begin() : leftBound;
        for (; resIt != rightBound; ++resIt) {
            assert((*resIt)->getLength() == (*resIt)->getSource()->getLength());
            oc = slowOverlap(resIt->get()->getTarget(), inputIt->get()->getTarget());
            // assert(oc == Same || oc == Disjoint || oc == AContainsB);
            if (oc == AContainsB) {
                clippedSegs.clear();
                clipAagainstB(*resIt, *inputIt, oc, clippedSegs);
                results.insert(clippedSegs.begin(), clippedSegs.end());
            }
        }
    }

    // 4) insert the input list
    results.insert(inputSegs.begin(), inputSegs.end());
}

// Map all segments from the input to any segments in the same genome
// that coalesce in or before the given "coalescence limit" genome.
// Destructive to any data in the input list.
static hal_size_t mapRecursiveParalogies(const Genome *srcGenome, list<MappedSegmentPtr> &input,
                                         list<MappedSegmentPtr> &results, const set<string> &namesOnPath,
                                         const Genome *coalescenceLimit, hal_size_t minLength) {
    if (input.empty()) {
        results = input;
        return 0;
    }

    const Genome *curGenome = (*input.begin())->getGenome();
    assert(curGenome != NULL);
    if (curGenome == coalescenceLimit) {
        results = input;
        return 0;
    }

    const Genome *nextGenome = curGenome->getParent();

    if (nextGenome == NULL) {
        throw hal_exception("Hit root genome when attempting to map paralogies");
    }
    list<MappedSegmentPtr> paralogs;
    // Map to any paralogs in the current genome.
    // FIXME: I think the original segments are included in this, which is a waste.
    list<MappedSegmentPtr>::iterator i = input.begin();
    for (; i != input.end(); ++i) {
        assert((*i)->getGenome() == curGenome);
        mapSelf(*i, paralogs, minLength);
    }

    if (nextGenome != coalescenceLimit) {
        list<MappedSegmentPtr> nextSegments;
        // Map all of the original segments (not the paralogs, which is a
        // waste) up to the next genome.
        i = input.begin();
        for (; i != input.end(); ++i) {
            assert((*i)->getGenome() == curGenome);
            mapUp(*i, nextSegments, true, minLength);
        }

        // Recurse on the mapped segments.
        mapRecursiveParalogies(srcGenome, nextSegments, results, namesOnPath, coalescenceLimit, minLength);
    }

    // Map all the paralogs we found in this genome back to the source.
    list<MappedSegmentPtr> paralogsMappedToSrc;
    mapRecursiveDown(paralogs, paralogsMappedToSrc, srcGenome, namesOnPath, false, minLength);

    results.splice(results.begin(), paralogsMappedToSrc);
    results.sort(MappedSegment::LessSourcePtr());
    results.unique(MappedSegment::EqualToPtr());
    return results.size();
}

static hal_size_t mapSource(const SegmentIterator *source, MappedSegmentSet &results, const Genome *tgtGenome,
                            const set<const Genome *> *genomesOnPath, bool doDupes, hal_size_t minLength,
                            const Genome *coalescenceLimit, const Genome *mrca) {
    assert(source != NULL);

    // FIXME: why does target start out as source??  This is all a bit clunky
    SegmentIteratorPtr startSourceSegIt;
    SegmentIteratorPtr startTargetSegIt;
    if (source->isTop()) {
        startSourceSegIt = dynamic_cast<const TopSegmentIterator *>(source)->clone();
        startTargetSegIt = dynamic_cast<const TopSegmentIterator *>(source)->clone();
    } else {
        startSourceSegIt = dynamic_cast<const BottomSegmentIterator *>(source)->clone();
        startTargetSegIt = dynamic_cast<const BottomSegmentIterator *>(source)->clone();
    }

    MappedSegmentPtr newMappedSeg(new MappedSegment(startSourceSegIt, startTargetSegIt));

    list<MappedSegmentPtr> input;
    input.push_back(newMappedSeg);
    list<MappedSegmentPtr> output;

    set<string> namesOnPath;
    assert(genomesOnPath != NULL);
    for (set<const Genome *>::const_iterator i = genomesOnPath->begin(); i != genomesOnPath->end(); ++i) {
        namesOnPath.insert((*i)->getName());
    }

    // FIXME: using multiple lists is probably much slower than just
    // reusing the results list over and over.
    list<MappedSegmentPtr> upResults;
    // Map all segments up to the MRCA of src and tgt.
    if (source->getGenome() != mrca) {
        mapRecursiveUp(input, upResults, mrca, minLength);
    } else {
        upResults = input;
    }

    list<MappedSegmentPtr> paralogResults;
    // Map to all paralogs that coalesce in or below the coalescenceLimit.
    if (mrca != coalescenceLimit && doDupes) {
        mapRecursiveParalogies(mrca, upResults, paralogResults, namesOnPath, coalescenceLimit, minLength);
    } else {
        paralogResults = upResults;
    }

    // Finally, map back down to the target genome.
    if (tgtGenome != mrca) {
        mapRecursiveDown(paralogResults, output, tgtGenome, namesOnPath, doDupes, minLength);
    } else {
        output = paralogResults;
    }

    list<MappedSegmentPtr>::iterator outIt = output.begin();
    for (; outIt != output.end(); ++outIt) {
        insertAndBreakOverlaps(*outIt, results);
    }

    return output.size();
}

hal_size_t hal::halMapSegment(const SegmentIterator *source, MappedSegmentSet &outSegments, const Genome *tgtGenome,
                              const set<const Genome *> *genomesOnPath, bool doDupes, hal_size_t minLength,
                              const Genome *coalescenceLimit, const Genome *mrca) {
    assert(tgtGenome != NULL);

    if (mrca == NULL) {
        set<const Genome *> inputSet;
        inputSet.insert(source->getGenome());
        inputSet.insert(tgtGenome);
        mrca = getLowestCommonAncestor(inputSet);
    }

    if (coalescenceLimit == NULL) {
        coalescenceLimit = mrca;
    }

    // Get the path from the coalescence limit to the target (necessary
    // for choosing which children to move through to get to the
    // target).
    set<const Genome *> pathSet;
    if (genomesOnPath == NULL) {
        set<const Genome *> inputSet;
        inputSet.insert(tgtGenome);
        inputSet.insert(mrca);
        getGenomesInSpanningTree(inputSet, pathSet);
        genomesOnPath = &pathSet;
    }

    hal_size_t numResults =
        mapSource(source, outSegments, tgtGenome, genomesOnPath, doDupes, minLength, coalescenceLimit, mrca);
    return numResults;
}

/* call main function with smart pointer */
hal_size_t hal::halMapSegmentSP(const SegmentIteratorPtr &source, MappedSegmentSet &outSegments, const Genome *tgtGenome,
                                const std::set<const Genome *> *genomesOnPath, bool doDupes, hal_size_t minLength,
                                const Genome *coalescenceLimit, const Genome *mrca) {
    return halMapSegment(source.get(), outSegments, tgtGenome, genomesOnPath, doDupes, minLength, coalescenceLimit, mrca);
}
