/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <limits>
#include <algorithm>
#include <limits>
#include <cassert>
#include "halMappedSegment.h"
#include "halTopSegmentIterator.h"
#include "halBottomSegmentIterator.h"

using namespace std;
using namespace hal;

MappedSegment::MappedSegment(
  SegmentIteratorPtr sourceSegIt,
  SegmentIteratorPtr targetSegIt)
  :
    _source(std::dynamic_pointer_cast<SegmentIterator>(sourceSegIt)),
    _target(std::dynamic_pointer_cast<SegmentIterator>(targetSegIt))
{
  assert(_source->getLength() == _target->getLength());
}

bool hal::MappedSegmentLess::operator()(const hal::MappedSegment& m1,
                                        const hal::MappedSegment& m2) const {
    return m1.lessThan(&m2);
}

bool hal::MappedSegmentLess::operator()(const hal::MappedSegmentPtr& m1,
                                        const hal::MappedSegmentPtr& m2) const {
    return m1->lessThan(m2.get());
}



//////////////////////////////////////////////////////////////////////////////
// MAPPED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
bool MappedSegment::lessThan(const MappedSegment* other) const
{
  assert(other != NULL);
  int res = fastComp(_target, other->_target);
  if (res == 0)
  {
    res = fastComp(_source, other->_source);
  }
  return res == -1;
}

bool MappedSegment::lessThanBySource(
  const MappedSegment* other) const
{
  assert(other != NULL);
  int res = fastComp(_source, other->_source);
  if (res == 0)
  {
    res = fastComp(_target, other->_target);
  }
  return res == -1;
}

bool MappedSegment::equals(const MappedSegment* other) const
{
  assert(other != NULL);
  int res = fastComp(_source, other->_source);
  if (res == 0)
  {
    res = fastComp(_target, other->_target);
  }
  return res == 0;
}

void MappedSegment::flip()
{
  swap(_source, _target);
}

void MappedSegment::fullReverse()
{
  _source->slice(_source->getEndOffset(), _source->getStartOffset());
  _source->toReverse();
  _target->slice(_target->getEndOffset(), _target->getStartOffset());
  _target->toReverse();
  assert(_source->getLength() == _target->getLength());
}

MappedSegment* MappedSegment::clone() const
{
    // FIXME: having both sourceCpySegIt and sourceSegIt seems pointless, same for target.
  SegmentIteratorPtr sourceCpySegIt;
  if (_source->isTop())
  {
    sourceCpySegIt = std::dynamic_pointer_cast<TopSegmentIterator>(_source)->clone();
  }
  else
  {
    sourceCpySegIt = std::dynamic_pointer_cast<BottomSegmentIterator>(_source)->clone();
  }
  SegmentIteratorPtr sourceSegIt = std::dynamic_pointer_cast<SegmentIterator>(sourceCpySegIt);

  SegmentIteratorPtr targetCopySegIt;
  if (_target->isTop())
  {
      targetCopySegIt = std::dynamic_pointer_cast<TopSegmentIterator>(_target)->clone();
  } else {
      targetCopySegIt = std::dynamic_pointer_cast<BottomSegmentIterator>(_target)->clone();
  }
  SegmentIteratorPtr targetSegIt = 
     std::static_pointer_cast<SegmentIterator>(targetCopySegIt);

  assert(sourceSegIt->getStartPosition() == _source->getStartPosition() &&
         sourceSegIt->getEndPosition() == _source->getEndPosition());
  assert(targetSegIt->getStartPosition() == _target->getStartPosition() &&
         targetSegIt->getEndPosition() == _target->getEndPosition());
  assert(_source->getLength() == _target->getLength());
  assert(sourceSegIt->getLength() == targetSegIt->getLength());

  MappedSegment* newSeg = new MappedSegment(sourceSegIt, targetSegIt);
  
  assert(newSeg->getStartPosition() == getStartPosition() &&
         newSeg->getEndPosition() == getEndPosition() &&
         newSeg->_source->getStartPosition() == _source->getStartPosition() &&
         newSeg->_source->getEndPosition() == _source->getEndPosition());
  assert(newSeg->_source.get() != _source.get() &&
         newSeg->_target.get() != _target.get());
  return newSeg;
}

bool MappedSegment::canMergeRightWith(
  const MappedSegmentPtr& nextSeg,
  const set<hal_index_t>* cutSet,
  const set<hal_index_t>* sourceCutSet) const
{
  //return false;
  bool ret = false;
  const SlicedSegment* ref = this->getSource();
  const SlicedSegment* nextRef = nextSeg->getSource();
#if 0 // FIXME: why is this here
  assert(ref->getReversed() == false);
  assert(nextRef->getReversed() == false);
#endif
  assert(ref->getSequence() == nextRef->getSequence());
  assert(this->getGenome() == nextSeg->getGenome());
  hal_index_t sourceCut;
  hal_index_t cut;

  if (this->getReversed() == nextSeg->getReversed() &&
      ref->getReversed() == nextRef->getReversed())
  {
    hal_index_t qdelta = NULL_INDEX;
    hal_index_t rdelta = NULL_INDEX;
    if (this->getReversed() == false && ref->getReversed() == false)
    {
      qdelta = nextSeg->getStartPosition() - this->getEndPosition();
      rdelta = nextRef->getStartPosition() - ref->getEndPosition();
      cut = this->getEndPosition();
      sourceCut = ref->getEndPosition();
    }
    else if (this->getReversed() == true && ref->getReversed() == true)
    {
      qdelta = nextSeg->getEndPosition() - this->getStartPosition();
      rdelta = nextRef->getEndPosition() - ref->getStartPosition();
      cut = this->getStartPosition();
      sourceCut = ref->getStartPosition();
    }
    else if (this->getReversed() == false && ref->getReversed() == true)
    {
      qdelta = nextSeg->getStartPosition() - this->getEndPosition();
      rdelta = ref->getEndPosition() - nextRef->getStartPosition();
      cut = this->getEndPosition();
      sourceCut = nextRef->getStartPosition();
    }
    else  
    {
      assert(this->getReversed() == true && ref->getReversed() == false);
      qdelta = nextSeg->getEndPosition() - this->getStartPosition();
      rdelta = ref->getStartPosition() - nextRef->getEndPosition();
      cut = this->getStartPosition();
      sourceCut = nextRef->getEndPosition();
    }
    //assert(qdelta >= 0);
    ret = qdelta == 1 && rdelta == 1;
  }
  
  if (ret == true)
  {
    if (sourceCutSet != NULL && 
        sourceCutSet->find(sourceCut) != sourceCutSet->end())
    {
      ret = false;
    }
    else if (cutSet != NULL && 
             cutSet->find(cut) != cutSet->end())
    {
      ret = false;
    }
  }
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

int MappedSegment::fastComp(const SegmentIteratorPtr& s1, 
                            const SegmentIteratorPtr& s2)
{
  // compare without accessing anything from disk (ie using only index
  // and offset)
  int res = 0;
  assert(s1->getGenome() == s2->getGenome());
  if (s1->isTop() != s2->isTop())
  {
    res = boundComp(s1, s2);
    if (res == 0)
    {
      res = slowComp(s1, s2);
    }
  }
  else
  {
    if (s1->getArrayIndex() < s2->getArrayIndex())
    {
      res = -1;
    }
    else if (s1->getArrayIndex() > s2->getArrayIndex())
    {
      res = 1;
    }
    else
    { 
      hal_offset_t so1 = s1->getStartOffset();
      hal_offset_t eo1 = s1->getEndOffset();
      if (s1->getReversed())
      {
        swap(so1, eo1);
      }
      hal_offset_t so2 = s2->getStartOffset();
      hal_offset_t eo2 = s2->getEndOffset();
      if (s2->getReversed())
      {
        swap(so2, eo2);
      }
      if (so1 < so2)
      {
        res = -1;
      }
      else if (so1 > so2)
      {
        res = 1;
      }
      else if (eo1 > eo2)
      {
        res = -1;
      }
      else if (eo1 < eo2)
      {
        res = 1;
      }
    }
  }
  assert (res == slowComp(s1, s2));
  return res;
}

int MappedSegment::boundComp(const SegmentIteratorPtr& s1, 
                             const SegmentIteratorPtr& s2)
{
  int res = 0;
  bool flip = s2->getReversed();
  if (flip)
  {
    s2->toReverse();
  }
  
  if (s1->isTop() && !s2->isTop())
  {
    BottomSegmentIteratorPtr bot =
        std::dynamic_pointer_cast<BottomSegmentIterator>(s2);
    hal_index_t lb = bot->bs()->getTopParseIndex();
    hal_index_t ub = lb;
    if ((hal_size_t)bot->getArrayIndex() <
        bot->getGenome()->getNumBottomSegments()-1)
    {
      bot = bot->clone();
      bot->slice(0,0);
      bot->toRight();
      ub = bot->bs()->getTopParseIndex();
    }
    if (s1->getArrayIndex() < lb)
    {
      res = -1;
    }
    else if (s1->getArrayIndex() > ub)
    {
      res = 1;
    }
  }
  else if (!s1->isTop() && s2->isTop())
  {
    TopSegmentIteratorPtr top =
        std::dynamic_pointer_cast<TopSegmentIterator>(s2);
    hal_index_t lb = top->ts()->getBottomParseIndex();
    hal_index_t ub = lb;
    if ((hal_size_t)top->getArrayIndex() < 
        top->getGenome()->getNumTopSegments()-1)
    {
      top = top->clone();
      top->slice(0,0);
      top->toRight();
      ub = top->ts()->getBottomParseIndex();
    }
    if (s1->getArrayIndex() < lb)
    {
      res = -1;
    }
    else if (s1->getArrayIndex() > ub)
    {
      res = 1;
    }
  }

  if (flip)
  {
    s2->toReverse();
  }

  return res;
}

int MappedSegment::slowComp(const SegmentIteratorPtr& s1, 
                            const SegmentIteratorPtr& s2)
{
  assert(s1->getGenome() == s2->getGenome());
  int res = 0;
  hal_index_t sp1 = s1->getStartPosition();
  hal_index_t ep1 = s1->getEndPosition();
  hal_index_t sp2 = s2->getStartPosition();
  hal_index_t ep2 = s2->getEndPosition();
  if (s1->getReversed())
  {
    swap(sp1, ep1);
  }
  if (s2->getReversed())
  {
    swap(sp2, ep2);
  }
  if (sp1 < sp2)
  {
    res = -1;
  }
  else if (sp1 > sp2)
  {
    res = 1;
  }
  else if (ep1 < ep2)
  {
    res = -1;
  }
  else if (ep1 > ep2)
  {
    res = 1;
  }
  return res;
}


//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void MappedSegment::setArrayIndex(Genome* genome, 
                                  hal_index_t arrayIndex)
{
  _target->setArrayIndex(genome, arrayIndex);
}

const Genome* MappedSegment::getGenome() const
{
  return _target->getGenome();
}

Genome* MappedSegment::getGenome()
{
  return const_cast<Genome*>(_target->getGenome());
}

const Sequence* MappedSegment::getSequence() const
{
  return _target->getSequence();
}

hal_index_t MappedSegment::getStartPosition() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getStartPosition();
}

hal_index_t MappedSegment::getEndPosition() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getEndPosition();
}

hal_size_t MappedSegment::getLength() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getLength();
}

void MappedSegment::getString(string& outString) const
{
  _target->getString(outString);
}

void MappedSegment::setCoordinates(hal_index_t startPos, 
                                          hal_size_t length)
{
  throw hal_exception("MappedSegment::setCoordinates not imeplemented");
}

hal_index_t MappedSegment::getArrayIndex() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getArrayIndex();
}

bool MappedSegment::leftOf(hal_index_t genomePos) const
{
  return _target->leftOf(genomePos);
}

bool MappedSegment::rightOf(hal_index_t genomePos) const
{
  return _target->rightOf(genomePos);
}

bool MappedSegment::overlaps(hal_index_t genomePos) const
{
  return _target->overlaps(genomePos);
}

bool MappedSegment::isFirst() const
{
  return _target->isFirst();
}

bool MappedSegment::isLast() const
{
  return _target->isLast();
}

bool MappedSegment::isMissingData(double nThreshold) const
{
  return _target->isMissingData(nThreshold);
}

bool MappedSegment::isTop() const
{
  return _target->isTop();
}

void MappedSegment::print(ostream& os) const
{
  os << "Mapped Segment:\n";
  os << "Source: ";
  if (_source.get() == NULL)
  {
    os << "NULL";
  }
  else
  {
    os << *_source;
  }
  os << "\nTarget: ";
  if (_target.get() == NULL)
  {
    os << "NULL";
  }
  else
  {
    os << *_target;
  }
}

//////////////////////////////////////////////////////////////////////////////
// SLICED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////

void MappedSegment::toReverse()
{
  assert(_target->getLength() == _source->getLength());
  _target->toReverse();
}

void MappedSegment::toReverseInPlace()
{
  assert(_target->getLength() == _source->getLength());
  _target->toReverseInPlace();
}

hal_offset_t MappedSegment::getStartOffset() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getStartOffset();
}

hal_offset_t MappedSegment::getEndOffset() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getEndOffset();
}

void MappedSegment::slice(hal_offset_t startOffset ,
                                 hal_offset_t endOffset )
{
  assert(_source->getLength() == _target->getLength());
  hal_index_t startDelta = startOffset - _target->getStartOffset();
  hal_index_t endDelta = endOffset - _target->getEndOffset();
  _target->slice(startOffset, endOffset);
  _source->slice(_source->getStartOffset() + startDelta,
                 _source->getEndOffset() + endDelta);
  assert(_source->getLength() == _target->getLength());
}

bool MappedSegment::getReversed() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getReversed();
}

