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

MafExport::MafExport() : _maxRefGap(0), _noDupes(false)
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
  writeHeader();

  ColumnIteratorConstPtr colIt = seq->getColumnIterator(&targets,
                                                        _maxRefGap, 
                                                        startPosition,
                                                        lastPosition,
                                                        _noDupes,
                                                        _noAncestors);
  _mafBlock.initBlock(colIt, _ucscNames);
  assert(_mafBlock.canAppendColumn(colIt) == true);
  _mafBlock.appendColumn(colIt);
  size_t numBlocks = 0;
  while (colIt->lastColumn() == false)
  {
    colIt->toRight();
    if (_mafBlock.canAppendColumn(colIt) == false)
    {
      // erase empty entries from the column.  helps when there are 
      // millions of sequences (ie from fastas with lots of scaffolds)
      if (numBlocks++ % 1000 == 0)
      {
        colIt->defragment();
      }

      mafStream << _mafBlock << '\n';
      _mafBlock.initBlock(colIt, _ucscNames);
      assert(_mafBlock.canAppendColumn(colIt) == true);
    }
    _mafBlock.appendColumn(colIt);
  }
  mafStream << _mafBlock << endl;
}

