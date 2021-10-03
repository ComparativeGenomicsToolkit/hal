/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafExport.h"
#include <cassert>
#include <deque>

using namespace std;
using namespace hal;

void MafExport::writeHeader() {
    assert(_mafStream != NULL);
    // sometimes tellp() returns -1
    if (_mafStream->tellp() <= streampos(0)) {
        *_mafStream << "##maf version=1 scoring=N/A\n"
                    << "# hal " << _alignment->getNewickTree() << endl
                    << endl;
    }
}

void MafExport::convertSequence(ostream &mafStream, AlignmentConstPtr alignment, const Sequence *seq, hal_index_t startPosition,
                                hal_size_t length, const set<const Genome *> &targets) {
    assert(seq != NULL);
    if (startPosition >= (hal_index_t)seq->getSequenceLength() ||
        (hal_size_t)startPosition + length > seq->getSequenceLength()) {
        throw hal_exception("Invalid range specified for convertGenome");
    }
    if (length == 0) {
        length = seq->getSequenceLength() - startPosition;
    }
    if (length == 0) {
        throw hal_exception("Cannot convert zero length sequence");
    }
    hal_index_t lastPosition = startPosition + (hal_index_t)(length - 1);

    _mafStream = &mafStream;
    _alignment = alignment;
    if (!_append) {
        writeHeader();
    }

    ColumnIteratorPtr colIt = seq->getColumnIterator(&targets, _maxRefGap, startPosition, lastPosition, _noDupes, _noAncestors,
                                                     false, // reverseStrand,
                                                     true,  // unique
                                                     _onlyOrthologs);

    hal_size_t appendCount = 0;
    if (_unique == false || colIt->isCanonicalOnRef() == true) {
        _mafBlock.initBlock(colIt, _ucscNames, _printTree);
        assert(_mafBlock.canAppendColumn(colIt) == true);
        _mafBlock.appendColumn(colIt);
        ++appendCount;
    }
    size_t numBlocks = 0;
    while (colIt->lastColumn() == false) {
        colIt->toRight();
        if (_unique == false || colIt->isCanonicalOnRef() == true) {
            if (appendCount == 0) {
                _mafBlock.initBlock(colIt, _ucscNames, _printTree);
                assert(_mafBlock.canAppendColumn(colIt) == true);
            }
            if (_mafBlock.canAppendColumn(colIt) == false) {
                // erase empty entries from the column.  helps when there are
                // millions of sequences (ie from fastas with lots of scaffolds)
                if (numBlocks++ % 1000 == 0) {
                    colIt->defragment();
                }
                if ((appendCount > 0) and (_keepEmptyRefBlocks or (not _mafBlock.referenceIsAllGaps()))) {
                    mafStream << _mafBlock << '\n';
                }
                _mafBlock.initBlock(colIt, _ucscNames, _printTree);
                assert(_mafBlock.canAppendColumn(colIt) == true);
            }
            _mafBlock.appendColumn(colIt);
            ++appendCount;
        }
    }
    // if nothing was ever added (seems to happen in corner case where
    // all columns violate unique), mafBlock ostream operator will crash
    // so we do following check
    if ((appendCount > 0) and (_keepEmptyRefBlocks or (not _mafBlock.referenceIsAllGaps()))) {
        mafStream << _mafBlock << endl;
    }
}

void MafExport::convertEntireAlignment(ostream &mafStream, AlignmentConstPtr alignment) {
    hal_size_t appendCount = 0;
    size_t numBlocks = 0;

    _mafStream = &mafStream;
    _alignment = alignment;

    writeHeader();

    // Load in all leaves from alignment
    vector<const Genome *> leafGenomes = getLeafGenomes(alignment.get());

    ColumnIterator::VisitCache visitCache;
    // Go through all the genomes one by one, and spit out any columns
    // they participate in that we haven't seen.
    for (hal_size_t i = 0; i < leafGenomes.size(); i++) {
        const Genome *genome = leafGenomes[i];
        ColumnIteratorPtr colIt = genome->getColumnIterator(NULL, 0, 0, NULL_INDEX, _noDupes, _noAncestors,
                                                            false, // reverseStrand
                                                            true,  // unique
                                                            _onlyOrthologs);
        colIt->setVisitCache(&visitCache);
        // So that we don't accidentally visit the first column if it's
        // already been visited.
        colIt->toSite(0, genome->getSequenceLength() - 1);
        for (;;) {
            if (appendCount == 0) {
                _mafBlock.initBlock(colIt, _ucscNames, _printTree);
                assert(_mafBlock.canAppendColumn(colIt) == true);
            }
            if (_mafBlock.canAppendColumn(colIt) == false) {
                // erase empty entries from the column.  helps when there are
                // millions of sequences (ie from fastas with lots of scaffolds)
                if (numBlocks++ % 1000 == 0) {
                    colIt->defragment();
                }
                if (appendCount > 0) {
                    mafStream << _mafBlock << '\n';
                }
                _mafBlock.initBlock(colIt, _ucscNames, _printTree);
                assert(_mafBlock.canAppendColumn(colIt) == true);
            }
            _mafBlock.appendColumn(colIt);
            appendCount++;

            if (colIt->lastColumn()) {
                // Have to break here because otherwise
                // colIt->toRight() will crash.
                break;
            }
            colIt->toRight();
        }
        // Copy over the updated visit cache information. This is a
        // deep copy, so it's slow, but necessary to preserve the
        // column iterator ownership of the visit cache
        visitCache.clear();
        ColumnIterator::VisitCache *newVisitCache = colIt->getVisitCache();
        for (ColumnIterator::VisitCache::iterator it = newVisitCache->begin(); it != newVisitCache->end(); it++) {
            visitCache[it->first] = new PositionCache(*it->second);
        }
    }

    // if nothing was ever added (seems to happen in corner case where
    // all columns violate unique), mafBlock ostream operator will crash
    // so we do following check
    if (appendCount > 0) {
        mafStream << _mafBlock << endl;
    }
}
