/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hal.h"
#include "defaultMappedSegment.h"
#include "defaultTopSegmentIterator.h"
#include "defaultBottomSegmentIterator.h"

using namespace std;
using namespace hal;

DefaultMappedSegment::DefaultMappedSegment(
  DefaultSegmentIteratorConstPtr source,
  DefaultSegmentIteratorConstPtr target)
  :
  _source(source),
  _target(target)
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
    new DefaultMappedSegment(
      startSource.downCast<DefaultSegmentIteratorConstPtr>(), 
      startTarget.downCast<DefaultSegmentIteratorConstPtr>()));
  
  vector<DefaultMappedSegmentConstPtr> input(1, newMappedSeg);
  vector<DefaultMappedSegmentConstPtr> output;
  mapRecursive(source, input, output, tgtGenome, genomesOnPath, 
               doDupes);
  results.insert(results.end(), output.begin(), output.end());
  return output.size();
}

hal_size_t DefaultMappedSegment::mapRecursive(
  const DefaultSegmentIterator* source,
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

  if (!inputPtr->empty())
  {
    srcGenome = inputPtr->at(0)->getSource()->getGenome();
    genome = inputPtr->at(0)->getGenome();
  }

  if (doDupes == true)
  {
    outputPtr->clear();
    for (size_t i = 0; i < inputPtr->size(); ++i)
    {
      mapSelf(inputPtr->at(i), *outputPtr);
    }
    inputPtr->clear();
    swap(inputPtr, outputPtr);
  }

  if (!inputPtr->empty())     
  {
    const Genome* parentGenome = genome->getParent();
    if (parentGenome != NULL &&
        parentGenome != srcGenome && (
          parentGenome == tgtGenome || 
          genomesOnPath->find(parentGenome) != genomesOnPath->end()))
    {
      outputPtr->clear();
      for (size_t i = 0; i < inputPtr->size(); ++i)
      {
        mapUp(inputPtr->at(i), *outputPtr);
      }
      inputPtr->clear();
      swap(inputPtr, outputPtr);
      mapRecursive(source, *inputPtr, *outputPtr, tgtGenome, genomesOnPath, 
                   doDupes);
    }
    else
    {
      for (hal_size_t child = 0; child < genome->getNumChildren(); ++child)
      {
        // note that this code is potentially unfriendly to the 
        // inmemory option where entire child genomes get loaded
        // when they may not be needed.  but since the column iterator
        // does the same thing, we don't worry for now
        const Genome* childGenome = genome->getChild(child);
        if (childGenome != srcGenome && (
              childGenome == tgtGenome ||
              genomesOnPath->find(childGenome) != genomesOnPath->end()))
        {
          outputPtr->clear();
          for (size_t i = 0; i < inputPtr->size(); ++i)
          {
            mapDown(inputPtr->at(i), child, *outputPtr);
          }
          inputPtr->clear();
          swap(inputPtr, outputPtr);
          mapRecursive(source, *inputPtr, *outputPtr, tgtGenome, genomesOnPath, 
                       doDupes);
          break;
        }
      }
    }
  }
  
  // could potentially save this copy but dont care for now
  if (outputPtr != &results)
  {
    results = *outputPtr;
  }
  return results.size();
}

hal_size_t DefaultMappedSegment::mapUp(DefaultMappedSegmentConstPtr mappedSeg, 
                                       vector<DefaultMappedSegmentConstPtr>& results)
{
  hal_size_t added = 0;
  if (mappedSeg->isTop() == true)
  {
    const Genome* parent = mappedSeg->getGenome()->getParent();
    assert(parent != NULL);
    BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
    SegmentIteratorConstPtr target = mappedSeg->_target;
    TopSegmentIteratorConstPtr top = 
       target.downCast<TopSegmentIteratorConstPtr>();
    bottom->toParent(top);
    results.push_back(mappedSeg);
    ++added;
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
        new DefaultMappedSegment(
          topNew.downCast<DefaultSegmentIteratorConstPtr>(), 
          bottomNew.downCast<DefaultSegmentIteratorConstPtr>()));
      added += mapUp(newMappedSeg, results);
      top->toRight(rightCutoff);
    } while (top->getEndPosition() != rightCutoff);
  }
  return added;
}

hal_size_t DefaultMappedSegment::mapDown(
  DefaultMappedSegmentConstPtr mappedSeg, 
  hal_size_t childIndex,
  vector<DefaultMappedSegmentConstPtr>& results)
{
  hal_size_t added = 0;
  if (mappedSeg->isTop() == false)
  {
    const Genome* child = mappedSeg->getGenome()->getChild(childIndex);
    assert(child != NULL);
    TopSegmentIteratorConstPtr top = child->getTopSegmentIterator();
    SegmentIteratorConstPtr target = mappedSeg->_target;
    BottomSegmentIteratorConstPtr bottom = 
       target.downCast<BottomSegmentIteratorConstPtr>();
    top->toChild(bottom, childIndex);
    results.push_back(mappedSeg);
    ++added;
  }
  else
  {
    hal_index_t rightCutoff = mappedSeg->getEndPosition();
    BottomSegmentIteratorConstPtr bottom = 
       mappedSeg->getGenome()->getBottomSegmentIterator();
    SegmentIteratorConstPtr target = mappedSeg->_target;
    TopSegmentIteratorConstPtr top = 
       target.downCast<TopSegmentIteratorConstPtr>();
    bottom->toParseDown(top);
    do
    {
      BottomSegmentIteratorConstPtr bottomNew = bottom->copy();
      TopSegmentIteratorConstPtr topNew = top->copy();
      topNew->toParseUp(bottom);
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(
          bottomNew.downCast<DefaultSegmentIteratorConstPtr>(), 
          topNew.downCast<DefaultSegmentIteratorConstPtr>()));
      added += mapUp(newMappedSeg, results);
      bottom->toRight(rightCutoff);
    } while (bottom->getEndPosition() != rightCutoff);
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
      
      DefaultMappedSegmentConstPtr newMappedSeg(
        new DefaultMappedSegment(
          topCopy.downCast<DefaultSegmentIteratorConstPtr>(), 
          mappedSeg->_source));
      results.push_back(mappedSeg);
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
        new DefaultMappedSegment(
          topNew.downCast<DefaultSegmentIteratorConstPtr>(), 
          bottomNew.downCast<DefaultSegmentIteratorConstPtr>()));
      added += mapUp(newMappedSeg, results);
      top->toRight(rightCutoff);
    } while (top->getEndPosition() != rightCutoff);
  }
  return added;
}

// SEGMENT INTERFACE
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

// SLICED SEGMENT INTERFACE 
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
