/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halBlockLiftover.h"
#include "halBlockMapper.h"
#include "halBottomSegmentIterator.h"
#include "halMappedSegment.h"
#include "halSegmentMapper.h"
#include "halTopSegmentIterator.h"
#include <cassert>
#include <deque>
#include <list>

using namespace std;
using namespace hal;

BlockLiftover::BlockLiftover() : Liftover(), _useBatchPath(false), _closeGenomes(false) {
}

BlockLiftover::~BlockLiftover() {
}

void BlockLiftover::visitBegin() {
    if (_srcGenome->getNumTopSegments() > 0) {
        _refSeg = _srcGenome->getTopSegmentIterator();
        _lastIndex = (hal_index_t)_srcGenome->getNumTopSegments();
    } else {
        _refSeg = _srcGenome->getBottomSegmentIterator();
        _lastIndex = (hal_index_t)_srcGenome->getNumBottomSegments();
    }

    set<const Genome *> inputSet;
    inputSet.insert(_srcGenome);
    inputSet.insert(_tgtGenome);
    _mrca = getLowestCommonAncestor(inputSet);
    if (_coalescenceLimit == NULL) {
        _coalescenceLimit = _mrca;
    }

    inputSet.clear();
    inputSet.insert(_coalescenceLimit);
    inputSet.insert(_tgtGenome);
    getGenomesInSpanningTree(inputSet, _downwardPath);

    // Decide whether to use the batch path (genome-closing optimization).
    // Only possible when coalescenceLimit == MRCA.
    _useBatchPath = (_coalescenceLimit == _mrca);
    if (_useBatchPath) {
        _mrcaName = _mrca->getName();
        _downwardPathNames.clear();
        for (set<const Genome *>::const_iterator i = _downwardPath.begin(); i != _downwardPath.end(); ++i) {
            _downwardPathNames.insert((*i)->getName());
        }
    }
}

void BlockLiftover::liftIntervalBatchMap(hal_index_t globalStart, hal_index_t globalEnd, bool flip) {
    // Collect all source segments into a list of MappedSegments.
    list<MappedSegmentPtr> sourceSegs;
    while (_refSeg->getArrayIndex() < _lastIndex && _refSeg->getStartPosition() <= globalEnd) {
        if (flip == true) {
            _refSeg->toReverseInPlace();
        }

        SegmentIteratorPtr srcClone;
        SegmentIteratorPtr tgtClone;
        if (_refSeg->isTop()) {
            srcClone = dynamic_cast<const TopSegmentIterator *>(_refSeg.get())->clone();
            tgtClone = dynamic_cast<const TopSegmentIterator *>(_refSeg.get())->clone();
        } else {
            srcClone = dynamic_cast<const BottomSegmentIterator *>(_refSeg.get())->clone();
            tgtClone = dynamic_cast<const BottomSegmentIterator *>(_refSeg.get())->clone();
        }
        sourceSegs.push_back(MappedSegmentPtr(new MappedSegment(srcClone, tgtClone)));

        if (flip == true) {
            _refSeg->toReverseInPlace();
        }
        _refSeg->toRight(globalEnd);
    }

    // Map all segments through the tree in one batched pass.
    halMapSegmentBatch(sourceSegs, _mappedSegments, _tgtGenome, _downwardPathNames, _traverseDupes, 0, _mrcaName,
                       _srcGenome);
}

void BlockLiftover::liftInterval(BedList &mappedBedLines) {
    _mappedSegments.clear();
    hal_index_t globalStart = _bedLine._start + _srcSequence->getStartPosition();
    hal_index_t globalEnd = _bedLine._end - 1 + _srcSequence->getStartPosition();
    bool flip = _bedLine._strand == '-';

    _refSeg->toSite(globalStart, false);
    hal_offset_t startOffset = globalStart - _refSeg->getStartPosition();
    hal_offset_t endOffset = 0;
    if (globalEnd <= _refSeg->getEndPosition()) {
        endOffset = _refSeg->getEndPosition() - globalEnd;
    }
    _refSeg->slice(startOffset, endOffset);

    assert(_refSeg->getStartPosition() == globalStart);
    assert(_refSeg->getEndPosition() <= globalEnd);

    if (_useBatchPath && _closeGenomes) {
        // Batch path: collect all segments and map through the tree in one
        // pass, closing intermediate genomes to save memory.
        liftIntervalBatchMap(globalStart, globalEnd, flip);
    } else {
        // Per-segment path: required when coalescenceLimit != MRCA because
        // mapRecursiveParalogies re-traverses intermediate genomes.
        while (_refSeg->getArrayIndex() < _lastIndex && _refSeg->getStartPosition() <= globalEnd) {
            if (flip == true) {
                _refSeg->toReverseInPlace();
            }
            halMapSegment(_refSeg.get(), _mappedSegments, _tgtGenome, &_downwardPath, _traverseDupes, 0, _coalescenceLimit,
                          _mrca);
            if (flip == true) {
                _refSeg->toReverseInPlace();
            }
            _refSeg->toRight(globalEnd);
        }
    }

    vector<MappedSegmentPtr> fragments;
    MappedSegmentSet emptySet;
    set<hal_index_t> queryCutSet;
    set<hal_index_t> targetCutSet;

    for (MappedSegmentSet::iterator i = _mappedSegments.begin(); i != _mappedSegments.end(); ++i) {
        BlockMapper::extractSegment(i, emptySet, fragments, &_mappedSegments, targetCutSet, queryCutSet);

        const Sequence *seq = (*i)->getSequence();
        hal_size_t seqStart = seq->getStartPosition();
        mappedBedLines.push_back(_bedLine);
        BedLine &outBedLine = mappedBedLines.back();
        outBedLine._blocks.clear();
        outBedLine._chrName = seq->getName();
        outBedLine._start = min(min(fragments.front()->getStartPosition(), fragments.front()->getEndPosition()),
                                min(fragments.back()->getStartPosition(), fragments.back()->getEndPosition()));
        outBedLine._start -= seqStart;
        outBedLine._end = 1 + max(max(fragments.front()->getStartPosition(), fragments.front()->getEndPosition()),
                                  max(fragments.back()->getStartPosition(), fragments.back()->getEndPosition()));
        outBedLine._end -= seqStart;
        outBedLine._strand = (*i)->getReversed() ? '-' : '+';

        const SlicedSegment *srcFront = fragments.front()->getSource();
        const SlicedSegment *srcBack = fragments.back()->getSource();
        outBedLine._srcStart = min(min(srcFront->getStartPosition(), srcFront->getEndPosition()),
                                   min(srcBack->getStartPosition(), srcBack->getEndPosition()));
        outBedLine._srcStrand = srcFront->getReversed() ? '-' : '+';

        if (_bedLine._strand == '.') {
            outBedLine._strand = '.';
            outBedLine._srcStrand = '.';
        }

        assert(outBedLine._start < outBedLine._end);

        if (_outPSL == true && !fragments.empty()) {
            readPSLInfo(fragments, outBedLine);
        }
    }
}

void BlockLiftover::readPSLInfo(vector<MappedSegmentPtr> &fragments, BedLine &outBedLine) {
    const Sequence *srcSequence = fragments[0]->getSource()->getSequence();
    const Sequence *tSequence = fragments[0]->getSequence();

    outBedLine._psl.resize(1);
    PSLInfo &psl = outBedLine._psl[0];
    psl._matches = 0;
    psl._misMatches = 0;
    psl._repMatches = 0;
    psl._nCount = 0;
    psl._qNumInsert = 0;
    psl._qBaseInsert = 0;
    psl._tNumInsert = 0;
    psl._tBaseInsert = 0;
    psl._qSeqName = srcSequence->getName();
    psl._qSeqSize = srcSequence->getSequenceLength();
    psl._qStrand = fragments[0]->getSource()->getReversed() ? '-' : '+';
    assert(outBedLine._srcStart >= srcSequence->getStartPosition());
    psl._qChromOffset = srcSequence->getStartPosition();
    psl._qEnd = outBedLine._srcStart + (outBedLine._end - outBedLine._start);
    psl._tSeqSize = tSequence->getSequenceLength();
    psl._qBlockStarts.clear();

    string sBuf;
    string tBuf;
    for (size_t i = 0; i < fragments.size(); ++i) {
        assert(fragments[i]->getSource()->getReversed() == fragments[0]->getSource()->getReversed());
        assert(fragments[i]->getSource()->getSequence() == srcSequence);
        assert(fragments[i]->getSequence() == tSequence);

        fragments[i]->getSource()->getString(sBuf);
        fragments[i]->getString(tBuf);

        for (size_t j = 0; j < sBuf.length(); ++j) {
            if (sBuf[j] == tBuf[j]) {
                if (!isMasked(sBuf[j]) && !isMasked(tBuf[j])) {
                    ++psl._matches;
                } else {
                    ++psl._repMatches;
                }
            } else if (isMissingData(tBuf[j])) {
                ++psl._nCount;
            } else {
                ++psl._misMatches;
            }
        }
    }
}
