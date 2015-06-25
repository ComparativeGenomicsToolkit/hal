/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halMafExport.h"

using namespace std;
using namespace hal;

MafExport::MafExport() : _maxRefGap(0), _noDupes(false), _printTree(false)
{

}

MafExport::~MafExport()
{

}

void MafExport::setMaxRefGap(hal_size_t maxRefGap)
{
  _maxRefGap = maxRefGap;
}

void MafExport::setNoDupes(bool noDupes)
{
  _noDupes = noDupes;
}

void MafExport::setNoAncestors(bool noAncestors)
{
  _noAncestors = noAncestors;
}

void MafExport::setUcscNames(bool ucscNames)
{
  _ucscNames = ucscNames;
}

void MafExport::setUnique(bool unique)
{
  _unique = unique;
}

void MafExport::setAppend(bool append)
{
  _append = append;
}

void MafExport::setMaxBlockLength(hal_index_t maxLength)
{
  _mafBlock.setMaxLength(maxLength);
}

void MafExport::setPrintTree(bool printTree)
{
  _printTree = printTree;
}

void MafExport::setOnlyOrthologs(bool onlyOrthologs)
{
  _onlyOrthologs = onlyOrthologs;
}

void MafExport::writeHeader()
{
  assert(_mafStream != NULL);
  if (_mafStream->tellp() == streampos(0))
  {
    *_mafStream << "##maf version=1 scoring=N/A\n"
                << "# hal " << _alignment->getNewickTree() << endl << endl;
  }
}

void MafExport::convertSegmentedSequence(ostream& mafStream,
                                         AlignmentConstPtr alignment,
                                         const SegmentedSequence* seq,
                                         hal_index_t startPosition,
                                         hal_size_t length,
                                         const set<const Genome*>& targets)
{
  assert(seq != NULL);
  if (startPosition >= (hal_index_t)seq->getSequenceLength() ||
      (hal_size_t)startPosition + length > seq->getSequenceLength())
  {
    throw hal_exception("Invalid range specified for convertGenome");
  }
  if (length == 0)
  {
    length = seq->getSequenceLength() - startPosition;
  }
  if (length == 0)
  {
    throw hal_exception("Cannot convert zero length sequence");
  }
  hal_index_t lastPosition = startPosition + (hal_index_t)(length - 1);

  _mafStream = &mafStream;
  _alignment = alignment;
  if (!_append)
  {
    writeHeader();
  }

  ColumnIteratorConstPtr colIt = seq->getColumnIterator(&targets,
                                                        _maxRefGap, 
                                                        startPosition,
                                                        lastPosition,
                                                        _noDupes,
                                                        _noAncestors,
                                                        false, // reverseStrand,
                                                        true,  // unique
                                                        _onlyOrthologs);


  hal_size_t appendCount = 0;
  if (_unique == false || colIt->isCanonicalOnRef() == true)
  {
    _mafBlock.initBlock(colIt, _ucscNames, _printTree);
    assert(_mafBlock.canAppendColumn(colIt) == true);
    _mafBlock.appendColumn(colIt);
    ++appendCount;
  }
  size_t numBlocks = 0;
  while (colIt->lastColumn() == false)
  {
    colIt->toRight();
    if (_unique == false || colIt->isCanonicalOnRef() == true)
    {
      if (appendCount == 0)
      {
        _mafBlock.initBlock(colIt, _ucscNames, _printTree);
        assert(_mafBlock.canAppendColumn(colIt) == true);
      }
      if (_mafBlock.canAppendColumn(colIt) == false)
      {
        // erase empty entries from the column.  helps when there are 
        // millions of sequences (ie from fastas with lots of scaffolds)
        if (numBlocks++ % 1000 == 0)
        {
          colIt->defragment();
        }
        if (appendCount > 0)
        {
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
  if (appendCount > 0)
  {
    mafStream << _mafBlock << endl;
  }
}

void MafExport::convertEntireAlignment(ostream& mafStream,
                                       AlignmentConstPtr alignment)
{
    hal_size_t appendCount = 0;
    size_t numBlocks = 0;

    _mafStream = &mafStream;
    _alignment = alignment;

    writeHeader();

    // Load in all leaves from alignment
    vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
    vector<const Genome *> leafGenomes;
    for (hal_size_t i = 0; i < leafNames.size(); i++) {
        const Genome *genome = alignment->openGenome(leafNames[i]);
        assert(genome != NULL);
        leafGenomes.push_back(genome);
    }
    ColumnIterator::VisitCache visitCache;
    // Go through all the genomes one by one, and spit out any columns
    // they participate in that we haven't seen.
    for (hal_size_t i = 0; i < leafGenomes.size(); i++) {
        const Genome *genome = leafGenomes[i];
        ColumnIteratorConstPtr colIt = genome->getColumnIterator(NULL,
                                                                 0,
                                                                 0,
                                                                 NULL_INDEX,
                                                                 _noDupes,
                                                                 _noAncestors,
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
            if (_mafBlock.canAppendColumn(colIt) == false)
            {
                // erase empty entries from the column.  helps when there are 
                // millions of sequences (ie from fastas with lots of scaffolds)
                if (numBlocks++ % 1000 == 0)
                {
                    colIt->defragment();
                }
                if (appendCount > 0)
                {
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
        for(ColumnIterator::VisitCache::iterator it = newVisitCache->begin();
            it != newVisitCache->end(); it++) {
            visitCache[it->first] = new PositionCache(*it->second);
        }
    }

    // if nothing was ever added (seems to happen in corner case where
    // all columns violate unique), mafBlock ostream operator will crash
    // so we do following check
    if (appendCount > 0)
    {
        mafStream << _mafBlock << endl;
    }
}
