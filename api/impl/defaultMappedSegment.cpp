/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <limits>
#include "hal.h"
#include "defaultMappedSegment.h"
#include "defaultTopSegmentIterator.h"
#include "defaultBottomSegmentIterator.h"

using namespace std;
using namespace hal;

DefaultMappedSegment::DefaultMappedSegment(
  SegmentIteratorConstPtr source,
  SegmentIteratorConstPtr target)
  :
  _source(source.downCast<DefaultSegmentIteratorConstPtr>()),
  _target(target.downCast<DefaultSegmentIteratorConstPtr>())
{
  assert(_source->getLength() == _target->getLength());
}

DefaultMappedSegment::~DefaultMappedSegment()
{

}


//////////////////////////////////////////////////////////////////////////////
// MAPPED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
SlicedSegmentConstPtr DefaultMappedSegment::getSource() const
{
  return _source;
}

bool DefaultMappedSegment::lessThan(const MappedSegmentConstPtr& other) const
{
  DefaultMappedSegmentConstPtr od = const_cast<MappedSegmentConstPtr&>(
    other).downCast<DefaultMappedSegmentConstPtr>();
  assert(od.get() != NULL);
  int res = fastComp(_target, od->_target);
  if (res == 0)
  {
    res = fastComp(_source, od->_source);
  }
  return res == -1;
}

bool DefaultMappedSegment::lessThanBySource(
  const MappedSegmentConstPtr& other) const
{
  DefaultMappedSegmentConstPtr od = const_cast<MappedSegmentConstPtr&>(
    other).downCast<DefaultMappedSegmentConstPtr>();
  assert(od.get() != NULL);
  int res = fastComp(_source, od->_source);
  if (res == 0)
  {
    res = fastComp(_target, od->_target);
  }
  return res == -1;
}

bool DefaultMappedSegment::equals(const MappedSegmentConstPtr& other) const
{
  DefaultMappedSegmentConstPtr od = const_cast<MappedSegmentConstPtr&>(
    other).downCast<DefaultMappedSegmentConstPtr>();
  assert(od.get() != NULL);
  int res = fastComp(_source, od->_source);
  if (res == 0)
  {
    res = fastComp(_target, od->_target);
  }
  return res == 0;
}

void DefaultMappedSegment::flip() const
{
  swap(_source, _target);
}

void DefaultMappedSegment::fullReverse() const
{
  _source->slice(_source->getEndOffset(), _source->getStartOffset());
  _source->toReverse();
  _target->slice(_target->getEndOffset(), _target->getStartOffset());
  _target->toReverse();
  assert(_source->getLength() == _target->getLength());
}

MappedSegmentConstPtr DefaultMappedSegment::copy() const
{
  SegmentIteratorConstPtr srcCpy;
  if (_source->isTop())
  {
    srcCpy = _source.downCast<TopSegmentIteratorConstPtr>()->copy();
  }
  else
  {
    srcCpy = _source.downCast<BottomSegmentIteratorConstPtr>()->copy();
  }
  DefaultSegmentIteratorConstPtr srcIt = 
     srcCpy.downCast<DefaultSegmentIteratorConstPtr>();

  SegmentIteratorConstPtr tgtCpy;
  if (_target->isTop())
  {
    tgtCpy = _target.downCast<TopSegmentIteratorConstPtr>()->copy();
  }
  else
  {
    tgtCpy = _target.downCast<BottomSegmentIteratorConstPtr>()->copy();
  }
  DefaultSegmentIteratorConstPtr tgtIt = 
     tgtCpy.downCast<DefaultSegmentIteratorConstPtr>();

  assert(srcIt->getStartPosition() == _source->getStartPosition() &&
         srcIt->getEndPosition() == _source->getEndPosition());
  assert(tgtIt->getStartPosition() == _target->getStartPosition() &&
         tgtIt->getEndPosition() == _target->getEndPosition());
  assert(_source->getLength() == _target->getLength());
  assert(srcIt->getLength() == tgtIt->getLength());

  DefaultMappedSegment* newSeg = new DefaultMappedSegment(srcIt, tgtIt);
  
  assert(newSeg->getStartPosition() == getStartPosition() &&
         newSeg->getEndPosition() == getEndPosition() &&
         newSeg->_source->getStartPosition() == _source->getStartPosition() &&
         newSeg->_source->getEndPosition() == _source->getEndPosition());
  assert(newSeg->_source.get() != _source.get() &&
         newSeg->_target.get() != _target.get());
  return MappedSegmentConstPtr(newSeg);
}

bool DefaultMappedSegment::canMergeRightWith(
  const MappedSegmentConstPtr& next,
  const set<hal_index_t>* cutSet,
  const set<hal_index_t>* sourceCutSet) const
{
  //return false;
  bool ret = false;
  SlicedSegmentConstPtr ref = this->getSource();
  SlicedSegmentConstPtr nextRef = next->getSource();
//  assert(ref->getReversed() == false);
//  assert(nextRef->getReversed() == false);
  assert(ref->getSequence() == nextRef->getSequence());
  assert(this->getGenome() == next->getGenome());
  hal_index_t sourceCut;
  hal_index_t cut;

  if (this->getReversed() == next->getReversed() &&
      ref->getReversed() == nextRef->getReversed())
  {
    hal_index_t qdelta = NULL_INDEX;
    hal_index_t rdelta = NULL_INDEX;
    if (this->getReversed() == false && ref->getReversed() == false)
    {
      qdelta = next->getStartPosition() - this->getEndPosition();
      rdelta = nextRef->getStartPosition() - ref->getEndPosition();
      cut = this->getEndPosition();
      sourceCut = ref->getEndPosition();
    }
    else if (this->getReversed() == true && ref->getReversed() == true)
    {
      qdelta = next->getEndPosition() - this->getStartPosition();
      rdelta = nextRef->getEndPosition() - ref->getStartPosition();
      cut = this->getStartPosition();
      sourceCut = ref->getStartPosition();
    }
    else if (this->getReversed() == false && ref->getReversed() == true)
    {
      qdelta = next->getStartPosition() - this->getEndPosition();
      rdelta = ref->getEndPosition() - nextRef->getStartPosition();
      cut = this->getEndPosition();
      sourceCut = nextRef->getStartPosition();
    }
    else  
    {
      assert(this->getReversed() == true && ref->getReversed() == false);
      qdelta = next->getEndPosition() - this->getStartPosition();
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

int DefaultMappedSegment::fastComp(const DefaultSegmentIteratorConstPtr& s1, 
                                   const DefaultSegmentIteratorConstPtr& s2)
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

int DefaultMappedSegment::boundComp(const DefaultSegmentIteratorConstPtr& s1, 
                                    const DefaultSegmentIteratorConstPtr& s2)
{
  int res = 0;
  bool flip = s2->getReversed();
  if (flip)
  {
    s2->toReverse();
  }
  
  if (s1->isTop() && !s2->isTop())
  {
    BottomSegmentIteratorConstPtr bot = 
       const_cast<DefaultSegmentIteratorConstPtr&>(
         s2).downCast<BottomSegmentIteratorConstPtr>();
    hal_index_t lb = bot->getTopParseIndex();
    hal_index_t ub = lb;
    if ((hal_size_t)bot->getArrayIndex() <
        bot->getGenome()->getNumBottomSegments()-1)
    {
      bot = bot->copy();
      bot->slice(0,0);
      bot->toRight();
      ub = bot->getTopParseIndex();
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
    TopSegmentIteratorConstPtr top = 
       const_cast<DefaultSegmentIteratorConstPtr&>(
         s2).downCast<TopSegmentIteratorConstPtr>();
    hal_index_t lb = top->getBottomParseIndex();
    hal_index_t ub = lb;
    if ((hal_size_t)top->getArrayIndex() < 
        top->getGenome()->getNumTopSegments()-1)
    {
      top = top->copy();
      top->slice(0,0);
      top->toRight();
      ub = top->getBottomParseIndex();
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

int DefaultMappedSegment::slowComp(const DefaultSegmentIteratorConstPtr& s1, 
                                   const DefaultSegmentIteratorConstPtr& s2)
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

hal_size_t DefaultMappedSegment::map(const DefaultSegmentIterator* source,
                                     set<MappedSegmentConstPtr>& results,
                                     const Genome* tgtGenome,
                                     const set<const Genome*>* genomesOnPath,
                                     bool doDupes,
                                     hal_size_t minLength,
                                     const Genome *coalescenceLimit,
                                     const Genome *mrca)
{
  assert(source != NULL);
 
  SegmentIteratorConstPtr startSource;
  SegmentIteratorConstPtr startTarget;
  if (source->isTop())
  {
    startSource = 
       dynamic_cast<const DefaultTopSegmentIterator*>(source)->copy();
    startTarget = 
       dynamic_cast<const DefaultTopSegmentIterator*>(source)->copy();
  }
  else
  {
    startSource = 
       dynamic_cast<const DefaultBottomSegmentIterator*>(source)->copy();
    startTarget = 
       dynamic_cast<const DefaultBottomSegmentIterator*>(source)->copy();
  }
  
  DefaultMappedSegmentConstPtr newMappedSeg(
    new DefaultMappedSegment(startSource, startTarget)); 
  
  list<DefaultMappedSegmentConstPtr> input;
  input.push_back(newMappedSeg);
  list<DefaultMappedSegmentConstPtr> output;

  set<string> namesOnPath;
  assert(genomesOnPath != NULL);
  for (set<const Genome*>::const_iterator i = genomesOnPath->begin();
       i != genomesOnPath->end(); ++i)
  {
    namesOnPath.insert((*i)->getName());
  }

  // FIXME: using multiple lists is probably much slower than just
  // reusing the results list over and over.
  list<DefaultMappedSegmentConstPtr> upResults;
  // Map all segments up to the MRCA of src and tgt.
  if (source->getGenome() != mrca)
  {
    mapRecursiveUp(input, upResults, mrca, minLength);
  } else {
    upResults = input;
  }

  list<DefaultMappedSegmentConstPtr> paralogResults;
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

  list<DefaultMappedSegmentConstPtr>::iterator outIt = output.begin();
  for (; outIt != output.end(); ++outIt)
  {
    insertAndBreakOverlaps(*outIt, results);
  }

  return output.size();
}

// Map all segments from the input to any segments in the same genome
// that coalesce in or before the given "coalescence limit" genome.
// Destructive to any data in the input list.
hal_size_t DefaultMappedSegment::mapRecursiveParalogies(
  const Genome *srcGenome,
  list<DefaultMappedSegmentConstPtr>& input,
  list<DefaultMappedSegmentConstPtr>& results,
  const set<string>& namesOnPath,
  const Genome* coalescenceLimit,
  hal_size_t minLength)
{
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
  list<DefaultMappedSegmentConstPtr> paralogs;
  // Map to any paralogs in the current genome.
  // FIXME: I think the original segments are included in this, which is a waste.
  list<DefaultMappedSegmentConstPtr>::iterator i = input.begin();
  for (; i != input.end(); ++i)
  {
    assert((*i)->getGenome() == curGenome);
    mapSelf((*i), paralogs, minLength);
  }

  if (nextGenome != coalescenceLimit) {
    list<DefaultMappedSegmentConstPtr> nextSegments;
    // Map all of the original segments (not the paralogs, which is a
    // waste) up to the next genome.
    i = input.begin();
    for (; i != input.end(); ++i)
    {
      assert((*i)->getGenome() == curGenome);
      mapUp((*i), nextSegments, true, minLength);
    }

    // Recurse on the mapped segments.
    mapRecursiveParalogies(srcGenome, nextSegments, results, namesOnPath, coalescenceLimit, minLength);
  }

  // Map all the paralogs we found in this genome back to the source.
  list<DefaultMappedSegmentConstPtr> paralogsMappedToSrc;
  mapRecursiveDown(paralogs, paralogsMappedToSrc, srcGenome, namesOnPath, false, minLength);

  results.splice(results.begin(), paralogsMappedToSrc);
  results.sort(DefaultMappedSegment::LessSource());
  results.unique(DefaultMappedSegment::EqualTo());
  return results.size();
}

// Map the input segments up until reaching the target genome. If the
// target genome is below the source genome, fail miserably.
// Destructive to any data in the input or results list.
hal_size_t DefaultMappedSegment::mapRecursiveUp(
  list<DefaultMappedSegmentConstPtr>& input,
  list<DefaultMappedSegmentConstPtr>& results,
  const Genome* tgtGenome,
  hal_size_t minLength)
{
  list<DefaultMappedSegmentConstPtr>* inputPtr = &input;
  list<DefaultMappedSegmentConstPtr>* outputPtr = &results;

  if (inputPtr->empty() || (*inputPtr->begin())->getGenome() == tgtGenome)
  {
    results = *inputPtr;
    return 0;
  }

  const Genome *curGenome = (*inputPtr->begin())->getGenome();
  assert(curGenome != NULL);
  const Genome *nextGenome = curGenome->getParent();

  if (nextGenome == NULL)
  {
    stringstream ss;
    ss << "Reached top of tree when attempting to recursively map up from "
       << curGenome->getName() << " to " << tgtGenome->getName();
    throw hal_exception(ss.str());
  }
  
  // Map all segments to the parent.
  list<DefaultMappedSegmentConstPtr>::iterator i = inputPtr->begin();
  for (; i != inputPtr->end(); ++i)
  {
    assert((*i)->getGenome() == curGenome);
    mapUp((*i), *outputPtr, true, minLength);
  }
  
  if (nextGenome != tgtGenome)
  {
    // Continue the recursion.
    swap(inputPtr, outputPtr);
    outputPtr->clear();
    mapRecursiveUp(*inputPtr, *outputPtr, tgtGenome, minLength);
  }

  if (outputPtr != &results)
  {
    results = *outputPtr;
  }

  results.sort(DefaultMappedSegment::LessSource());
  results.unique(DefaultMappedSegment::EqualTo());
  return results.size();
}

// Map the input segments down until reaching the target genome. If the
// target genome is above the source genome, fail miserably.
// Destructive to any data in the input or results list.
hal_size_t DefaultMappedSegment::mapRecursiveDown(
  list<DefaultMappedSegmentConstPtr>& input,
  list<DefaultMappedSegmentConstPtr>& results,
  const Genome* tgtGenome,
  const set<string>& namesOnPath,
  bool doDupes,
  hal_size_t minLength)
{
  list<DefaultMappedSegmentConstPtr>* inputPtr = &input;
  list<DefaultMappedSegmentConstPtr>* outputPtr = &results;

  if (inputPtr->empty())
  {
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
  for (hal_size_t child = 0; 
       nextGenome == NULL && child < childNames.size(); ++child)
  {
    if (childNames[child] == tgtGenome->getName() ||
        namesOnPath.find(childNames[child]) != namesOnPath.end())
    {
      const Genome* childGenome = curGenome->getChild(child);
      nextGenome = childGenome;
      nextChildIndex = child;
    }
  }

  if (nextGenome == NULL)
  {
    stringstream ss;
    ss << "Could not find correct child that leads from "
       << curGenome->getName() << " to " << tgtGenome->getName();
    throw hal_exception(ss.str());
  }

  assert(nextGenome->getParent() == curGenome);

  // Map the actual segments down.
  list<DefaultMappedSegmentConstPtr>::iterator i = inputPtr->begin();
  for (; i != inputPtr->end(); ++i)
  {
    assert((*i)->getGenome() == curGenome);
    mapDown((*i), nextChildIndex, *outputPtr, minLength);
  }

  // Find paralogs.
  if (doDupes == true)
  {
    swap(inputPtr, outputPtr);
    outputPtr->clear();
    list<DefaultMappedSegmentConstPtr>::iterator i = inputPtr->begin();
    for (; i != inputPtr->end(); ++i)
    {
      assert((*i)->getGenome() == nextGenome);
      mapSelf((*i), *outputPtr, minLength);
    }
  }

  if (nextGenome != tgtGenome)
  {
    // Continue the recursion.
    swap(inputPtr, outputPtr);
    outputPtr->clear();
    mapRecursiveDown(*inputPtr, *outputPtr, tgtGenome, namesOnPath, doDupes, minLength);
  }

  if (outputPtr != &results)
  {
    results = *outputPtr;
  }

  results.sort(DefaultMappedSegment::LessSource());
  results.unique(DefaultMappedSegment::EqualTo());
  return results.size();
}

hal_size_t DefaultMappedSegment::mapUp(
  DefaultMappedSegmentConstPtr mappedSeg, 
  list<DefaultMappedSegmentConstPtr>& results,
  bool doDupes,
  hal_size_t minLength)
{
  const Genome* parent = mappedSeg->getGenome()->getParent();
  assert(parent != NULL);
  hal_size_t added = 0;
  if (mappedSeg->isTop() == true)
  {
    BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
    TopSegmentIteratorConstPtr top = mappedSeg->targetAsTop();
    if (top->hasParent() == true && top->getLength() >= minLength &&
        (doDupes == true || top->isCanonicalParalog() == true))
    {
      bottom->toParent(top);
      mappedSeg->_target = bottom.downCast<DefaultSegmentIteratorConstPtr>();
      results.push_back(mappedSeg);
      ++added;
    }
  }
  else
  {
    hal_index_t rightCutoff = mappedSeg->getEndPosition();
    BottomSegmentIteratorConstPtr bottom = mappedSeg->targetAsBottom();
    hal_index_t startOffset = (hal_index_t)bottom->getStartOffset();
    hal_index_t endOffset = (hal_index_t)bottom->getEndOffset();
    TopSegmentIteratorConstPtr top = 
       mappedSeg->getGenome()->getTopSegmentIterator();
    top->toParseUp(bottom);
    do
    {
      TopSegmentIteratorConstPtr topNew = top->copy();
      
      // we map the new target back to see how the offsets have 
      // changed.  these changes are then applied to the source segment
      // as deltas
      BottomSegmentIteratorConstPtr bottomBack = bottom->copy();
      bottomBack->toParseDown(topNew);
      hal_index_t startBack = (hal_index_t)bottomBack->getStartOffset();
      hal_index_t endBack = (hal_index_t)bottomBack->getEndOffset();
      assert(startBack >= startOffset);
      assert(endBack >= endOffset);
      SegmentIteratorConstPtr newSource = mappedSeg->sourceCopy();
      hal_index_t startDelta = startBack - startOffset;
      hal_index_t endDelta = endBack - endOffset;
      assert((hal_index_t)newSource->getLength() > startDelta + endDelta);
      newSource->slice(newSource->getStartOffset() + startDelta, 
                       newSource->getEndOffset() + endDelta);

      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(newSource, topNew));

      assert(newMappedSeg->isTop() == true);
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());

      added += mapUp(newMappedSeg, results, doDupes, minLength);
      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (top->getEndPosition() != rightCutoff)
      {
        top->toRight(rightCutoff);
      }
      else
      {
        break;
      }
    } 
    while (true);
  }
  return added;
}

hal_size_t DefaultMappedSegment::mapDown(
  DefaultMappedSegmentConstPtr mappedSeg, 
  hal_size_t childIndex,
  list<DefaultMappedSegmentConstPtr>& results,
  hal_size_t minLength)
{
  const Genome* child = mappedSeg->getGenome()->getChild(childIndex);
  assert(child != NULL);
  hal_size_t added = 0;
  if (mappedSeg->isTop() == false)
  {
    TopSegmentIteratorConstPtr top = child->getTopSegmentIterator();
    SegmentIteratorConstPtr target = mappedSeg->_target;
    BottomSegmentIteratorConstPtr bottom = 
       target.downCast<BottomSegmentIteratorConstPtr>();

    if (bottom->hasChild(childIndex) == true && bottom->getLength() >= minLength)
    {
      top->toChild(bottom, childIndex);
      mappedSeg->_target = top.downCast<DefaultSegmentIteratorConstPtr>();
      results.push_back(mappedSeg);
      ++added;
    }
  }
  else
  {
    hal_index_t rightCutoff = mappedSeg->getEndPosition();
    TopSegmentIteratorConstPtr top = mappedSeg->targetAsTop();
    hal_index_t startOffset = (hal_index_t)top->getStartOffset();
    hal_index_t endOffset = (hal_index_t)top->getEndOffset();
    BottomSegmentIteratorConstPtr bottom = 
       mappedSeg->getGenome()->getBottomSegmentIterator();
    bottom->toParseDown(top);
    do
    {
      BottomSegmentIteratorConstPtr bottomNew = bottom->copy();
      
      // we map the new target back to see how the offsets have 
      // changed.  these changes are then applied to the source segment
      // as deltas
      TopSegmentIteratorConstPtr topBack = top->copy();
      topBack->toParseUp(bottomNew);
      hal_index_t startBack = (hal_index_t)topBack->getStartOffset();
      hal_index_t endBack = (hal_index_t)topBack->getEndOffset();
      assert(startBack >= startOffset);
      assert(endBack >= endOffset);
      SegmentIteratorConstPtr newSource = mappedSeg->sourceCopy();
      hal_index_t startDelta = startBack - startOffset;
      hal_index_t endDelta = endBack - endOffset;
      assert((hal_index_t)newSource->getLength() > startDelta + endDelta);
      newSource->slice(newSource->getStartOffset() + startDelta, 
                       newSource->getEndOffset() + endDelta);
      
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(newSource, bottomNew));

      assert(newMappedSeg->isTop() == false);
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());

      added += mapDown(newMappedSeg, childIndex, results, minLength);

      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (bottom->getEndPosition() != rightCutoff)
      {
        bottom->toRight(rightCutoff);
      }
      else
      {
        break;
      }
    } 
    while (true);
  }
  return added;
}

hal_size_t DefaultMappedSegment::mapSelf(
  DefaultMappedSegmentConstPtr mappedSeg, 
  list<DefaultMappedSegmentConstPtr>& results,
  hal_size_t minLength)
{
  hal_size_t added = 0;
  if (mappedSeg->isTop() == true)
  {
    SegmentIteratorConstPtr target = mappedSeg->_target;
    SegmentIteratorConstPtr source = mappedSeg->_source;
    TopSegmentIteratorConstPtr top = 
       target.downCast<TopSegmentIteratorConstPtr>();
    TopSegmentIteratorConstPtr topCopy = top->copy();
    do
    {
      SegmentIteratorConstPtr newSource = source->isTop() ? 
         (SegmentIteratorConstPtr)source.downCast<
           TopSegmentIteratorConstPtr>()->copy() :
         (SegmentIteratorConstPtr)source.downCast<
           BottomSegmentIteratorConstPtr>()->copy();
      TopSegmentIteratorConstPtr newTop = topCopy->copy();
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(newSource, newTop));
      assert(newMappedSeg->getGenome() == mappedSeg->getGenome());
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());
      results.push_back(newMappedSeg);
      ++added;
      if (topCopy->hasNextParalogy())
      {
        topCopy->toNextParalogy();
      }
    } 
    while (topCopy->hasNextParalogy() == true && 
           topCopy->getLength() >= minLength &&
           topCopy->getArrayIndex() != top->getArrayIndex());
  }
  else if (mappedSeg->getGenome()->getParent() != NULL)
  {
    hal_index_t rightCutoff = mappedSeg->getEndPosition();
    BottomSegmentIteratorConstPtr bottom = mappedSeg->targetAsBottom();
    hal_index_t startOffset = (hal_index_t)bottom->getStartOffset();
    hal_index_t endOffset = (hal_index_t)bottom->getEndOffset();
    TopSegmentIteratorConstPtr top = 
       mappedSeg->getGenome()->getTopSegmentIterator();
    top->toParseUp(bottom);
    do
    {
      TopSegmentIteratorConstPtr topNew = top->copy();
      
      // we map the new target back to see how the offsets have 
      // changed.  these changes are then applied to the source segment
      // as deltas
      BottomSegmentIteratorConstPtr bottomBack = bottom->copy();
      bottomBack->toParseDown(topNew);
      hal_index_t startBack = (hal_index_t)bottomBack->getStartOffset();
      hal_index_t endBack = (hal_index_t)bottomBack->getEndOffset();
      assert(startBack >= startOffset);
      assert(endBack >= endOffset);
      SegmentIteratorConstPtr newSource = mappedSeg->sourceCopy();
      hal_index_t startDelta = startBack - startOffset;
      hal_index_t endDelta = endBack - endOffset;
      assert((hal_index_t)newSource->getLength() > startDelta + endDelta);
      newSource->slice(newSource->getStartOffset() + startDelta, 
                       newSource->getEndOffset() + endDelta);

      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(newSource, topNew));

      assert(newMappedSeg->isTop() == true);
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());

      added += mapSelf(newMappedSeg, results, minLength);
      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (top->getEndPosition() != rightCutoff)
      {
        top->toRight(rightCutoff);
      }
      else
      {
        break;
      }
    } 
    while (true);
  }
  return added;
}

DefaultMappedSegment::OverlapCat DefaultMappedSegment::slowOverlap(
  const SlicedSegmentConstPtr& sA, 
  const SlicedSegmentConstPtr& sB)
{
  hal_index_t startA = sA->getStartPosition();
  hal_index_t endA = sA->getEndPosition();
  hal_index_t startB = sB->getStartPosition();
  hal_index_t endB = sB->getEndPosition();
  if (startA > endA)
  {
    swap(startA, endA);
  }
  if (startB > endB)
  {
    swap(startB, endB);
  }

  if (endA < startB || startA > endB)
  {
    return Disjoint;
  }
  else if (startA == startB && endA == endB)
  {
    return Same;
  }
  else if (startA >= startB && endA <= endB)
  {
    return BContainsA;
  }
  else if (startB >= startA && endB <= endA)
  {
    return AContainsB;
  }
  else if (startA <= startB && endA < endB)
  {
    return AOverlapsLeftOfB;
  }
  assert (startB <= startA && endB < endA);
  return BOverlapsLeftOfA;
}

void DefaultMappedSegment::getOverlapBounds(
  const DefaultMappedSegmentConstPtr& seg, 
  set<MappedSegmentConstPtr>& results, 
  set<MappedSegmentConstPtr>::iterator& leftBound, 
  set<MappedSegmentConstPtr>::iterator& rightBound)
{
  if (results.size() <= 2)
  {
    leftBound = results.begin();
    rightBound = results.end();
  }
  else
  {
    set<MappedSegmentConstPtr>::iterator i = results.lower_bound(seg);
    leftBound = i;
    if (leftBound != results.begin())
    {
      --leftBound;
    }
    set<MappedSegmentConstPtr>::iterator iprev;
    set<MappedSegmentConstPtr>::key_compare resLess = results.key_comp();
    while (leftBound != results.begin())
    {
      iprev = leftBound;
      --iprev;
      if (leftBound == results.end() || !resLess(*iprev, *leftBound))
      {
        leftBound = iprev;
      }
      else
      {
        break;
      }
    }
    for (; leftBound != results.begin(); --leftBound)
    {
      if (leftBound != results.end() && 
          slowOverlap(seg, *leftBound) == Disjoint)
      {
        break;
      }
    }
    rightBound = i;
    if (rightBound != results.end())
    {
      for (++rightBound; rightBound != results.end(); ++rightBound)
      {
        if (slowOverlap(seg, *rightBound) == Disjoint)
        {
          break;
        }
      }
    }
  }
}

void 
DefaultMappedSegment::clipAagainstB(MappedSegmentConstPtr segA,
                                    MappedSegmentConstPtr segB,
                                    OverlapCat overlapCat,
                                    vector<MappedSegmentConstPtr>& clippedSegs)
{
  assert(overlapCat != Same && overlapCat != Disjoint && 
         overlapCat != BContainsA);

  hal_index_t startA = segA->getStartPosition();
  hal_index_t endA = segA->getEndPosition();
  hal_index_t startB = segB->getStartPosition();
  hal_index_t endB = segB->getEndPosition();
  if (startA > endA)
  {
    swap(startA, endA);
  }
  if (startB > endB)
  {
    swap(startB, endB);
  }
  MappedSegmentPtr left = segA;
  MappedSegmentPtr middle = segA->copy();
  MappedSegmentPtr right;

  hal_index_t startO = segA->getStartOffset();
  hal_index_t endO = segA->getEndOffset();
  hal_index_t length = segA->getLength();
  hal_index_t leftSize = std::max((hal_index_t)0, startB - startA);
  hal_index_t rightSize = std::max((hal_index_t)0, endA - endB);
  hal_index_t middleSize = length - leftSize - rightSize;
  if (rightSize > 0)
  {
    right = segA->copy();
  }

  assert (overlapCat == AOverlapsLeftOfB || overlapCat == BOverlapsLeftOfA ||
          overlapCat == AContainsB);
 
  assert(leftSize >= 0 && rightSize >=0 && middleSize >= 0);
  hal_index_t leftSlice = 0;
  hal_index_t rightSlice = 0;
  if (leftSize > 0)
  {
    leftSlice = 0;
    rightSlice = (length - leftSize);
    if (left->getReversed())
    {
      swap(leftSlice, rightSlice);
    }
    left->slice(startO + leftSlice, endO + rightSlice);
    assert(left->getLength() == (hal_size_t)leftSize);
    assert(min(left->getStartPosition(), left->getEndPosition()) == startA);
    assert(max(left->getStartPosition(), left->getEndPosition()) == startB - 1);
  }
  else
  {
    middle = segA;
  }
    
  leftSlice = leftSize;
  rightSlice = rightSize;
  if (middle->getReversed())
  {
    swap(leftSlice, rightSlice);
  }
  middle->slice(startO + leftSlice, endO + rightSlice);
  assert(middle->getLength() == (hal_size_t)middleSize);
  assert(min(middle->getStartPosition(), middle->getEndPosition()) == 
         max(startB, startA));
  assert(max(middle->getStartPosition(), middle->getEndPosition()) == 
         min(endB, endA));  
  if (middle.get() != segA.get())
  {
    assert(leftSize > 0);
    assert(middle->getLength() == middle->getSource()->getLength());
    clippedSegs.push_back(middle);
  }
  
  if (rightSize > 0)
  {
    leftSlice = leftSize + middleSize;
    rightSlice = 0;
    if (right->getReversed())
    {
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

void DefaultMappedSegment::insertAndBreakOverlaps(
  DefaultMappedSegmentConstPtr seg,
  set<MappedSegmentConstPtr>& results)
{
  // MappedSegmentConstPtr blin = seg->copy();
  // blin->slice(seg->getStartOffset(), 
  //             seg->getEndOffset());
  // results.insert(blin); return;
  assert(seg->getLength() == seg->getSource()->getLength());
  list<MappedSegmentConstPtr> inputSegs;
  vector<MappedSegmentConstPtr> clippedSegs;
  
  // 1) compute invariant range in set of candidate overalaps
  set<MappedSegmentConstPtr>::iterator leftBound;
  set<MappedSegmentConstPtr>::iterator rightBound;
  getOverlapBounds(seg, results, leftBound, rightBound);
  bool leftBegin = leftBound == results.begin();

  // 2) cut seg by each segment in range
  list<MappedSegmentConstPtr>::iterator inputIt;
  OverlapCat oc;  
  inputSegs.push_back(seg);
  set<MappedSegmentConstPtr>::iterator resIt;
  for (resIt = leftBound; resIt != rightBound; ++resIt)
  {
    for (inputIt = inputSegs.begin(); inputIt != inputSegs.end(); ++inputIt)
    {
      oc = slowOverlap(*inputIt, *resIt);
      if (oc == AContainsB || oc == AOverlapsLeftOfB || oc == BOverlapsLeftOfA)
      {
        clippedSegs.clear();
        clipAagainstB(*inputIt, *resIt, oc, clippedSegs);
        inputSegs.insert(inputSegs.end(), clippedSegs.begin(), 
                         clippedSegs.end());
      }
    }
  }
  
  // 3) cut results by input list
  for (inputIt = inputSegs.begin(); inputIt != inputSegs.end(); ++inputIt)
  {
    assert((*inputIt)->getLength() == (*inputIt)->getSource()->getLength());
    resIt = leftBegin ? results.begin() : leftBound;
    for (; resIt != rightBound; ++resIt)
    {
      assert((*resIt)->getLength() == (*resIt)->getSource()->getLength());
      oc = slowOverlap(*resIt, *inputIt);
      //assert(oc == Same || oc == Disjoint || oc == AContainsB);
      if (oc == AContainsB)
      {
        clippedSegs.clear();
        clipAagainstB(*resIt, *inputIt, oc, clippedSegs);
        results.insert(clippedSegs.begin(), clippedSegs.end());
      }
    }
  } 
  
  // 4) insert the input list
  results.insert(inputSegs.begin(), inputSegs.end());
}

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultMappedSegment::setArrayIndex(Genome* genome, 
                                         hal_index_t arrayIndex)
{
  _target->setArrayIndex(genome, arrayIndex);
}

void DefaultMappedSegment::setArrayIndex(const Genome* genome, 
                                         hal_index_t arrayIndex) const
{
  _target->setArrayIndex(genome, arrayIndex);
}

const Genome* DefaultMappedSegment::getGenome() const
{
  return _target->getGenome();
}

Genome* DefaultMappedSegment::getGenome()
{
  return const_cast<Genome*>(_target->getGenome());
}

const Sequence* DefaultMappedSegment::getSequence() const
{
  return _target->getSequence();
}

Sequence* DefaultMappedSegment::getSequence()
{
  return const_cast<Sequence*>(_target->getSequence());
}

hal_index_t DefaultMappedSegment::getStartPosition() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getStartPosition();
}

hal_index_t DefaultMappedSegment::getEndPosition() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getEndPosition();
}

hal_size_t DefaultMappedSegment::getLength() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getLength();
}

void DefaultMappedSegment::getString(string& outString) const
{
  _target->getString(outString);
}

void DefaultMappedSegment::setCoordinates(hal_index_t startPos, 
                                          hal_size_t length)
{
  throw hal_exception("DefaultMappedSegment::setCoordinates not imeplemented");
}

hal_index_t DefaultMappedSegment::getArrayIndex() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getArrayIndex();
}

bool DefaultMappedSegment::leftOf(hal_index_t genomePos) const
{
  return _target->leftOf(genomePos);
}

bool DefaultMappedSegment::rightOf(hal_index_t genomePos) const
{
  return _target->rightOf(genomePos);
}

bool DefaultMappedSegment::overlaps(hal_index_t genomePos) const
{
  return _target->overlaps(genomePos);
}

bool DefaultMappedSegment::isFirst() const
{
  return _target->isFirst();
}

bool DefaultMappedSegment::isLast() const
{
  return _target->isLast();
}

bool DefaultMappedSegment::isMissingData(double nThreshold) const
{
  return _target->isMissingData(nThreshold);
}

bool DefaultMappedSegment::isTop() const
{
  return _target->isTop();
}

hal_size_t DefaultMappedSegment::getMappedSegments(
  set<MappedSegmentConstPtr>& outSegments,
  const Genome* tgtGenome,
  const set<const Genome*>* genomesOnPath,
  bool doDupes,
  hal_size_t minLength,
  const Genome *coalescenceLimit,
  const Genome *mrca) const
{
  return _target->getMappedSegments(outSegments, tgtGenome, genomesOnPath,
                                    doDupes, minLength, coalescenceLimit,
                                    mrca);
}

void DefaultMappedSegment::print(ostream& os) const
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

void DefaultMappedSegment::toReverse() const
{
  assert(_target->getLength() == _source->getLength());
  _target->toReverse();
}

void DefaultMappedSegment::toReverseInPlace() const
{
  assert(_target->getLength() == _source->getLength());
  _target->toReverseInPlace();
}

hal_offset_t DefaultMappedSegment::getStartOffset() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getStartOffset();
}

hal_offset_t DefaultMappedSegment::getEndOffset() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getEndOffset();
}

void DefaultMappedSegment::slice(hal_offset_t startOffset ,
                                 hal_offset_t endOffset ) const
{
  assert(_source->getLength() == _target->getLength());
  hal_index_t startDelta = startOffset - _target->getStartOffset();
  hal_index_t endDelta = endOffset - _target->getEndOffset();
  _target->slice(startOffset, endOffset);
  _source->slice(_source->getStartOffset() + startDelta,
                 _source->getEndOffset() + endDelta);
  assert(_source->getLength() == _target->getLength());
}

bool DefaultMappedSegment::getReversed() const
{
  assert(_target->getLength() == _source->getLength());
  return _target->getReversed();
}
