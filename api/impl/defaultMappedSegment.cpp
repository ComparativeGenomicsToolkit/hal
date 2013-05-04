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

//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

hal_size_t DefaultMappedSegment::map(const DefaultSegmentIterator* source,
                                     vector<MappedSegmentConstPtr>& results,
                                     const Genome* tgtGenome,
                                     const set<const Genome*>* genomesOnPath,
                                     bool doDupes)
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
  
  vector<DefaultMappedSegmentConstPtr> input(1, newMappedSeg);
  vector<DefaultMappedSegmentConstPtr> output;
  mapRecursive(NULL, input, output, tgtGenome, genomesOnPath, doDupes);
  results.insert(results.end(), output.begin(), output.end());
  return output.size();
}

hal_size_t DefaultMappedSegment::mapRecursive(
  const Genome* prevGenome,
  vector<DefaultMappedSegmentConstPtr>& input,
  vector<DefaultMappedSegmentConstPtr>& results,
  const Genome* tgtGenome,
  const set<const Genome*>* genomesOnPath,
  bool doDupes)
{
  vector<DefaultMappedSegmentConstPtr>* inputPtr = &input;
  vector<DefaultMappedSegmentConstPtr>* outputPtr = &results;
  
  const Genome* srcGenome = NULL;
  const Genome* genome = NULL;
  const Genome* nextGenome = NULL;
  hal_size_t nextChildIndex = numeric_limits<hal_size_t>::max();

  if (!inputPtr->empty())
  {
    srcGenome = inputPtr->at(0)->getSource()->getGenome();
    genome = inputPtr->at(0)->getGenome();

    const Genome* parentGenome = genome->getParent();
    if (parentGenome != NULL &&
        parentGenome != prevGenome && (
          parentGenome == tgtGenome || 
          genomesOnPath->find(parentGenome) != genomesOnPath->end()))
    {
      nextGenome = parentGenome;
    }
    for (hal_size_t child = 0; 
         nextGenome == NULL && child < genome->getNumChildren(); ++child)
    {
      // note that this code is potentially unfriendly to the 
      // inmemory option where entire child genomes get loaded
      // when they may not be needed.  but since the column iterator
      // does the same thing, we don't worry for now
      const Genome* childGenome = genome->getChild(child);
      if (childGenome != prevGenome && (
            childGenome == tgtGenome ||
            genomesOnPath->find(childGenome) != genomesOnPath->end()))
      {
        nextGenome = childGenome;
        nextChildIndex = child;
      }
    }
  }
  
  if (nextGenome != NULL)
  {
    outputPtr->clear();

    if (nextGenome == genome->getParent())
    {
      for (size_t i = 0; i < inputPtr->size(); ++i)
      {
        mapUp(inputPtr->at(i), *outputPtr);
      }
    }
    else
    {
      for (size_t i = 0; i < inputPtr->size(); ++i)
      {
        mapDown(inputPtr->at(i), nextChildIndex, *outputPtr);
      }
      if (doDupes == true)
      {   
        size_t outputSize = outputPtr->size();
        for (size_t i = 0; i <outputSize; ++i)
        {
          mapSelf(outputPtr->at(i), *outputPtr);
        }
      }
    }
    swap(inputPtr, outputPtr);
    assert(genome != NULL);
    
    mapRecursive(genome, *inputPtr, *outputPtr, tgtGenome, genomesOnPath, 
                 doDupes);
  }
  else
  {
    swap(inputPtr, outputPtr);
  }
  
  // could potentially save this copy but dont care for now
  if (outputPtr != &results)
  {
    results = *outputPtr;
  }
  return results.size();
}

hal_size_t DefaultMappedSegment::mapUp(
  DefaultMappedSegmentConstPtr mappedSeg, 
  vector<DefaultMappedSegmentConstPtr>& results)
{
  const Genome* parent = mappedSeg->getGenome()->getParent();
  assert(parent != NULL);
  hal_size_t added = 0;
  if (mappedSeg->isTop() == true)
  {
    BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
    TopSegmentIteratorConstPtr top = mappedSeg->targetAsTop();
    if (top->hasParent() == true)
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
      assert(newSource->getLength() > startDelta + endDelta);
      newSource->slice(newSource->getStartOffset() + startDelta, 
                       newSource->getEndOffset() + endDelta);

      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(newSource, topNew));

      assert(newMappedSeg->isTop() == true);
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());

      added += mapUp(newMappedSeg, results);
      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (top->getEndPosition() != rightCutoff)
      {
        top->toRight(rightCutoff);
      }
    } 
    while (top->getEndPosition() != rightCutoff);
  }
  return added;
}

hal_size_t DefaultMappedSegment::mapDown(
  DefaultMappedSegmentConstPtr mappedSeg, 
  hal_size_t childIndex,
  vector<DefaultMappedSegmentConstPtr>& results)
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
    if (bottom->hasChild(childIndex) == true)
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

      added += mapDown(newMappedSeg, childIndex, results);

      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (bottom->getEndPosition() != rightCutoff)
      {
        bottom->toRight(rightCutoff);
      }
    } 
    while (bottom->getEndPosition() != rightCutoff);
  }
  return added;
}

hal_size_t DefaultMappedSegment::mapSelf(
  DefaultMappedSegmentConstPtr mappedSeg, 
  vector<DefaultMappedSegmentConstPtr>& results)
{
  hal_size_t added = 0;
  if (mappedSeg->isTop() == true)
  {
    SegmentIteratorConstPtr target = mappedSeg->_target;
    TopSegmentIteratorConstPtr top = 
       target.downCast<TopSegmentIteratorConstPtr>();
    TopSegmentIteratorConstPtr topCopy = top->copy();
    while (topCopy->hasNextParalogy() == true)
    {
      topCopy->toNextParalogy();
      if (topCopy->getArrayIndex() == top->getArrayIndex())
      {
        break;
      }
      SegmentIteratorConstPtr source = mappedSeg->_source;
      TopSegmentIteratorConstPtr newTop = topCopy->copy();
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(
          mappedSeg->_source,
          newTop.downCast<DefaultSegmentIteratorConstPtr>()));
      results.push_back(newMappedSeg);
      ++added;
    }
  }
  else
  {
    hal_index_t rightCutoff = mappedSeg->getEndPosition();
    TopSegmentIteratorConstPtr top = 
       mappedSeg->getGenome()->getTopSegmentIterator();
    SegmentIteratorConstPtr target = mappedSeg->_target;
    BottomSegmentIteratorConstPtr bottom = 
       target.downCast<BottomSegmentIteratorConstPtr>();
    top->toParseUp(bottom);
    do
    {
      TopSegmentIteratorConstPtr topNew = top->copy();
      BottomSegmentIteratorConstPtr bottomNew = bottom->copy();
      bottomNew->toParseDown(top);
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(bottomNew, topNew));
      assert(newMappedSeg->isTop() == true);
      assert(newMappedSeg->getSource()->getGenome() == 
             mappedSeg->getSource()->getGenome());
      added += mapUp(newMappedSeg, results);
      
      // stupid that we have to make this check but odn't want to 
      // make fundamental api change now
      if (top->getEndPosition() != rightCutoff)
      {
        top->toRight(rightCutoff);
      }
    } 
    while (top->getEndPosition() != rightCutoff);
  }
  return added;
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
  return _target->getStartPosition();
}

hal_index_t DefaultMappedSegment::getEndPosition() const
{
  return _target->getEndPosition();
}

hal_size_t DefaultMappedSegment::getLength() const
{
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
  vector<MappedSegmentConstPtr>& outSegments,
  const Genome* tgtGenome,
  const set<const Genome*>* genomesOnPath,
  bool doDupes) const
{
  return _target->getMappedSegments(outSegments, tgtGenome, genomesOnPath,
                                    doDupes);
}

//////////////////////////////////////////////////////////////////////////////
// SLICED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////

void DefaultMappedSegment::toReverse() const
{
  _target->toReverse();
}

hal_offset_t DefaultMappedSegment::getStartOffset() const
{
  return _target->getStartOffset();
}

hal_offset_t DefaultMappedSegment::getEndOffset() const
{
  return _target->getEndOffset();
}

void DefaultMappedSegment::slice(hal_offset_t startOffset ,
                                 hal_offset_t endOffset ) const
{
  throw hal_exception("DefaultMappedSegment::slice not implemented");
}

bool DefaultMappedSegment::getReversed() const
{
  return _target->getReversed();
}
