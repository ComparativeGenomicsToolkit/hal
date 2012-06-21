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

MafExport::MafExport()
{

}

MafExport::~MafExport()
{

}

void MafExport::writeHeader()
{
  assert(_mafStream != NULL);
  if (_mafStream->tellp() == 0)
  {
    *_mafStream << "##maf version=1 scoring=N/A\n"
                << "# hal " << _alignment->getNewickTree() << endl << endl;
  }
}

void MafExport::convertGenome(ostream& mafStream,
                              AlignmentConstPtr alignment,
                              const Genome* genome,
                              hal_index_t startPosition,
                              hal_size_t length,
                              const Genome* root)
{
  if (length == 0)
  {
    length = genome->getSequenceLength();
  }
  if ((hal_size_t)startPosition + length > genome->getSequenceLength())
  {
    throw hal_exception("Invalid range specified for convertGenome");
  }
  _mafStream = &mafStream;
  _alignment = alignment;
  writeHeader();
  hal_index_t pos = startPosition;
  hal_size_t doneLen = 0;

  while (doneLen < length)
  {
    const Sequence* sequence = genome->getSequenceBySite(pos);
    assert(sequence != NULL);
    hal_index_t seqStart = startPosition - sequence->getStartPosition();
    hal_size_t curLen = min(length - doneLen, 
                            sequence->getSequenceLength() - startPosition);
    assert(curLen > 0);
    convertSequence(mafStream, alignment, sequence, seqStart, curLen, root);
    doneLen += curLen;
  }
}

void MafExport::convertSequence(ostream& mafStream,
                                AlignmentConstPtr alignment,
                                const Sequence* sequence,
                                hal_index_t startPosition,
                                hal_size_t length,
                                const Genome* root)
{
  if (length == 0)
  {
    length = sequence->getSequenceLength();
  }
  if ((hal_size_t)startPosition + length > sequence->getSequenceLength())
  {
    throw hal_exception("Invalid range specified for convertSequence");
  }
  _mafStream = &mafStream;
  _alignment = alignment;
  writeHeader();

  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(root,
                                                             0, 
                                                             startPosition);
  ColumnIteratorConstPtr colEnd = 
     sequence->getColumnEndIterator(startPosition + length - 1);

  _mafBlock.setFirstColumn(colIt);
  assert(_mafBlock.canAppendColumn(colIt) == true);
  
  while (colIt < colEnd)
  {
    if (_mafBlock.canAppendColumn(colIt) == false)
    {
      mafStream << _mafBlock << endl << endl;
      _mafBlock.setFirstColumn(colIt);
      assert(_mafBlock.canAppendColumn(colIt) == true);
    }
    _mafBlock.appendColumn(colIt);
  }
  mafStream << _mafBlock << endl << endl;
}

