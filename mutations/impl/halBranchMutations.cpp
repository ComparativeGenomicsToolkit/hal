/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <locale>
#include "halBranchMutations.h"

using namespace std;
using namespace hal;

const string BranchMutations::inversionBedTag = "V";
const string BranchMutations::insertionBedTag = "I";
const string BranchMutations::deletionBedTag = "D";
const string BranchMutations::transpositionBedTag = "P";
const string BranchMutations::duplicationBedTag = "U";
const string BranchMutations::gapInsertionBedTag = "GI";
const string BranchMutations::gapDeletionBedTag = "GD";
string BranchMutations::substitutionBedTag(char parent, char child)
{
  assert(parent == toupper(parent));
  assert(child == toupper(child));
  assert(isNucleotide(parent) && parent != 'N' && 
         isNucleotide(child) && child != 'N');
  return string("S_") + parent + child;
}


BranchMutations::BranchMutations()
{

}

BranchMutations::~BranchMutations()
{

}

void BranchMutations::analyzeBranch(AlignmentConstPtr alignment,
                                    hal_size_t gapThreshold,
                                    double nThreshold,
                                    ostream* refBedStream,
                                    ostream* parentBedStream,
                                    ostream* snpBedStream,
                                    const Genome* reference,
                                    hal_index_t startPosition,
                                    hal_size_t length)
{
  assert(reference != NULL);
  if (length == 0)
  {
    throw hal_exception("Cannot convert zero length sequence");
  }
  if (reference->getParent() == NULL)
  {
    throw hal_exception("Reference genome must have parent");
  }

  _reference = reference;
  _start = startPosition;
  _length = length;
  _alignment = alignment;
  _maxGap = gapThreshold;
  _nThreshold = nThreshold;
  _refStream = refBedStream;
  _parentStream = parentBedStream;
  _snpStream = snpBedStream;
  _refName = _reference->getName();
  _parName = _reference->getParent()->getName();

  hal_index_t end = startPosition + (hal_index_t)length - 1;

  _top  = reference->getTopSegmentIterator();
  _top->toSite(startPosition);
  _bottom1 = reference->getParent()->getBottomSegmentIterator();
  _bottom2 = reference->getParent()->getBottomSegmentIterator();
  
  _rearrangement = reference->getRearrangement(_top->getArrayIndex(),
                                               _maxGap, _nThreshold);
  
  writeHeaders();

  do {
    _sequence = _reference->getSequenceBySite(
      _rearrangement->getLeftBreakpoint()->getStartPosition());
    switch (_rearrangement->getID())
    {
    case Rearrangement::Deletion:
      writeDeletion();
      break;
    case Rearrangement::Inversion:     
    case Rearrangement::Transposition:
    case Rearrangement::Insertion:
      writeInsertionOrInversion();
      break;
    case Rearrangement::Duplication:
      writeDuplication();
    default:
      break;
    }
    writeSubstitutions();
    writeGapInsertions();
  }
  while (_rearrangement->identifyNext() == true && 
         _rearrangement->getLeftBreakpoint()->getStartPosition() <= end);
  
  // kind of stupid, but we do a second pass to get the gapped deletions
  _top  = reference->getTopSegmentIterator();
  _top->toSite(startPosition);
  _rearrangement = reference->getRearrangement(_top->getArrayIndex(), 0,
                                               _nThreshold, true);   
  do {
    switch (_rearrangement->getID())
    {
    case Rearrangement::Deletion:
      writeGapDeletion();
    default:
      break;
    }
  }
  while (_rearrangement->identifyNext() == true && 
         _rearrangement->getLeftBreakpoint()->getStartPosition() <= end);
}

void BranchMutations::writeInsertionOrInversion()
{
  assert(_refStream != NULL);
  hal_size_t startPos = _rearrangement->getLeftBreakpoint()->getStartPosition();
  hal_size_t endPos = _rearrangement->getRightBreakpoint()->getEndPosition();
  
  assert(startPos < _start + _length && endPos >= _start);

  if (startPos < _start)
  {
    startPos = _start;
  }
  if (endPos >= _start + _length)
  {
    endPos = _start + _length - 1;
  }
  
  *_refStream << _sequence->getName() << '\t' 
              << startPos - _sequence->getStartPosition() << '\t'
              << endPos + 1 - _sequence->getStartPosition() << '\t';

  if (_rearrangement->getID() ==  Rearrangement::Inversion)
  {
    *_refStream << inversionBedTag << '\t';
  }
  else if (_rearrangement->getID() ==  Rearrangement::Insertion)
  {
    *_refStream << insertionBedTag << '\t';
  }
  else
  {
    assert(_rearrangement->getID() ==  Rearrangement::Transposition);
    *_refStream << transpositionBedTag << '\t';
  }
  
  *_refStream << _parName << '\t' << _refName << '\n';
}

void BranchMutations::writeSubstitutions()
{
  if (_snpStream == NULL)
  {
    return;
  }
  string tstring, bstring;
  hal_size_t pos;
  _top->copy(_rearrangement->getLeftBreakpoint());
  assert(_top->getReversed() == false);
  do {
    if (_top->hasParent())
    {
      _bottom1->toParent(_top);
      _top->getString(tstring);
      _bottom1->getString(bstring);
      assert(tstring.length() == bstring.length());

      for (hal_index_t i = 0; i < (hal_index_t)tstring.length(); ++i)
      {
        pos = i + _top->getStartPosition();
        char c = toupper(tstring[i]);
        char p = toupper(bstring[i]);

        if (pos >= _start && pos < _start + _length && 
            c != 'N' && p != 'N' && c != p)
        {        
          *_snpStream << _sequence->getName() << '\t'
                      << pos - _sequence->getStartPosition() << '\t' 
                      << pos + 1 - _sequence->getStartPosition() << '\t'
                      << substitutionBedTag(p, c) << '\t'
                      << _parName << '\t' << _refName<< '\n';
        }
      }
    }
    _top->toRight();
  } 
  while (_top->getArrayIndex() < 
         _rearrangement->getRightBreakpoint()->getArrayIndex());
}

void BranchMutations::writeGapInsertions()
{
  hal_size_t startPos, endPos;
  _top->copy(_rearrangement->getLeftBreakpoint());
  do {
    if (_top->hasParent() == false && _top->getLength() <= _maxGap)
    {
      startPos = _top->getStartPosition();
      endPos = _top->getEndPosition();
      assert(startPos < _start + _length && endPos >= _start);
      if (startPos < _start)
      {
        startPos = _start;
      }
      if (endPos >= _start + _length)
      {
        endPos = _start + _length - 1;
      }

      *_refStream << _sequence->getName() << '\t' 
                  << startPos - _sequence->getStartPosition() << '\t'
                  << endPos + 1 - _sequence->getStartPosition() << '\t'
                  << gapInsertionBedTag << '\t'
                  << _parName << '\t' << _refName << '\n';      
    }
    _top->toRight();
  }
  while (_top->getArrayIndex() < 
         _rearrangement->getRightBreakpoint()->getArrayIndex());
  
}

void BranchMutations::writeDeletion()
{
  if (_parentStream == NULL)
  {
    return;
  }

  pair<hal_index_t, hal_index_t> pos = _rearrangement->getDeletedRange();
  assert((pos.second - pos.first) + 1 > (hal_index_t)_maxGap);
  const Sequence* seq = _reference->getParent()->getSequenceBySite(pos.first);
  assert(seq != NULL);
  
  *_parentStream << seq->getName() << '\t' 
                 << pos.first - seq->getStartPosition() << '\t'
                 << pos.second + 1 - seq->getStartPosition() << '\t'
                 << deletionBedTag << '\t'
                 << _parName << '\t' << _refName << '\n'; 
}

void BranchMutations::writeGapDeletion()
{
  if (_parentStream == NULL)
  {
    return;
  }

  pair<hal_index_t, hal_index_t> pos = _rearrangement->getDeletedRange();
  if ((pos.second - pos.first) + 1 < (hal_index_t)_maxGap)
  {
    const Sequence* seq = _reference->getParent()->getSequenceBySite(pos.first);
    assert(seq != NULL);
  
    *_parentStream << seq->getName() << '\t' 
                   << pos.first - seq->getStartPosition() << '\t'
                   << pos.second + 1 - seq->getStartPosition() << '\t'
                   << gapDeletionBedTag << '\t'
                   << _parName << '\t' << _refName << '\n';
  }  
}

void BranchMutations::writeDuplication()
{
  if (_parentStream == NULL)
  {
    return;
  }
  pair<hal_index_t, hal_index_t> pos = _rearrangement->getDuplicatedRange();
  
  const Sequence* seq = _reference->getParent()->getSequenceBySite(pos.first);
  assert(seq != NULL);
  
  *_parentStream << seq->getName() << '\t' 
                 << pos.first - seq->getStartPosition() << '\t'
                 << pos.second + 1 - seq->getStartPosition() << '\t'
                 << duplicationBedTag << '\t'
                 << _parName << '\t' << _refName << '\n';
}

void BranchMutations::writeHeaders()
{
  string header("#Sequence\tStart\tEnd\tMutationID\tParentGenome\tChildGenome\n"
                "#I=Insertion D=Deletion GI(D)=GapInsertion(GapDeletion) "
                "V=Inversion P=Transposition U=Duplication\n");
  *_refStream << header;
  if (_parentStream)
  {
    *_parentStream << header;
  }
  if (_snpStream)
  {
    *_snpStream << header;
  }
  
}
