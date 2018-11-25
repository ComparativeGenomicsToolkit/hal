/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <deque>
#include <vector>
#include <iostream>
#include <cassert>
#include "halCommon.h"
#include "halValidate.h"
#include "halTopSegment.h"
#include "halTopSegmentIterator.h"
#include "halBottomSegment.h"
#include "halBottomSegmentIterator.h"
#include "halSequenceIterator.h"
#include "halGenome.h"
#include "halDnaIterator.h"

using namespace std;
using namespace hal;

// current implementation is poor and hacky.  should fix up to 
// use iterators to properly scan the segments. 

void hal::validateBottomSegment(const BottomSegment* bottomSegment)
{
  const Genome* genome = bottomSegment->getGenome();
  hal_index_t index = bottomSegment->getArrayIndex();
  if (index < 0 || index >= (hal_index_t)genome->getSequenceLength())
  {
      throw hal_exception("Bottom segment out of range " + std::to_string(index) + " in genome "
                        + genome->getName());
  }
  
  if (bottomSegment->getLength() < 1)
  {
      throw hal_exception("Bottom segment " + std::to_string(index)  + " in genome " + genome->getName()
                          + " has length 0 which is not currently supported");
  }

  hal_size_t numChildren = bottomSegment->getNumChildren();
  for (hal_size_t child = 0; child < numChildren; ++child)
  {
    const Genome* childGenome = genome->getChild(child);
    const hal_index_t childIndex = bottomSegment->getChildIndex(child);
    if (childGenome != NULL && childIndex != NULL_INDEX)
    {
      if (childIndex >= (hal_index_t)childGenome->getNumTopSegments())
      {
          throw hal_exception("Child " + std::to_string(child) + " index " + std::to_string(childIndex) + " of segment "
                            + std::to_string(bottomSegment->getArrayIndex()) + " out of range in genome "
                            + childGenome->getName());
      }
      TopSegmentIteratorPtr topSegmentIteratr = 
         childGenome->getTopSegmentIterator(childIndex);
      const TopSegment* childSegment = topSegmentIteratr->getTopSegment();
      if (childSegment->getLength() != bottomSegment->getLength())
      {
          throw hal_exception("Child " + std::to_string(child) + " with index " 
                            + std::to_string(childSegment->getArrayIndex())
                            + " and start position " + std::to_string(childSegment->getStartPosition())
                            + " and sequence " + childSegment->getSequence()->getName()
                            + " has length " + std::to_string(childSegment->getLength())
                            + " but parent with index " + std::to_string(bottomSegment->getArrayIndex() )
                            + " and start position " + std::to_string(bottomSegment->getStartPosition())
                            + " in sequence " + bottomSegment->getSequence()->getName() 
                            + " has length " + std::to_string(bottomSegment->getLength()));
      }
      if (childSegment->getNextParalogyIndex() == NULL_INDEX &&
          childSegment->getParentIndex() != bottomSegment->getArrayIndex())
      {
        throw hal_exception("Parent / child index mismatch:\n" 
                            + genome->getName() + "[" + std::to_string(bottomSegment->getArrayIndex()) + "]"
                            + " links to " + childGenome->getName() + "[" + std::to_string(childIndex)
                            + "] but \n"
                            + childGenome->getName() + "[" + std::to_string(childSegment->getArrayIndex())
                            + "] links to " + genome->getName() + "[" 
                            + std::to_string(childSegment->getParentIndex()) + "]");
      }
      if (childSegment->getParentReversed() != 
          bottomSegment->getChildReversed(child))
      {
          throw hal_exception("parent / child reversal mismatch (parent=" +
                              genome->getName() + " parentSegNum=" + std::to_string(bottomSegment->getArrayIndex()) + " child=" +
                              childGenome->getName() + " childSegNum=" + std::to_string(childSegment->getArrayIndex()) + ")");
      }
    }
  }

  const hal_index_t parseIndex = bottomSegment->getTopParseIndex();
  if (parseIndex == NULL_INDEX)
  {
    if (genome->getParent() != NULL)
    {
        throw hal_exception("Bottom segment " + std::to_string(bottomSegment->getArrayIndex()) + " in genome "
                          + genome->getName() + " has null parse index");
    }
  }
  else
  {
    if (parseIndex >= (hal_index_t)genome->getNumTopSegments())
    {
        throw hal_exception("BottomSegment " + std::to_string(bottomSegment->getArrayIndex()) + " in genome "
                            + genome->getName() + " has parse index " + std::to_string(parseIndex)
                            + " greater than the number of top segments, " 
                            + std::to_string(genome->getNumTopSegments()));
    }
    TopSegmentIteratorPtr parseIterator = 
       genome->getTopSegmentIterator(parseIndex);
    const TopSegment* parseSegment = parseIterator->getTopSegment();
    hal_offset_t parseOffset = bottomSegment->getTopParseOffset();
    if (parseOffset >= parseSegment->getLength())
    {
        throw hal_exception("BottomSegment " + std::to_string(bottomSegment->getArrayIndex()) + " in genome "
                            + genome->getName() + " has parse offset, " + std::to_string(parseOffset)
                            + ", greater than the length of the segment, "
                            + std::to_string(parseSegment->getLength()));
    }
    if ((hal_index_t)parseOffset + parseSegment->getStartPosition() != 
        bottomSegment->getStartPosition())
    {
      throw hal_exception("parse index broken in bottom segment in genome " +
                          genome->getName());
                          
    }
  }
}

void hal::validateTopSegment(const TopSegment* topSegment)
{
  const Genome* genome = topSegment->getGenome();
  hal_index_t index = topSegment->getArrayIndex();
  if (index < 0 || index >= (hal_index_t)genome->getSequenceLength())
  {
      throw hal_exception("Segment out of range " + std::to_string(index) + " in genome "
                        + genome->getName());
  }

  if (topSegment->getLength() < 1)
  {
      throw hal_exception("Top segment " + std::to_string(index)  + " in genome " + genome->getName()
                          + " has length 0 which is not currently supported");
  }

  const Genome* parentGenome = genome->getParent();
  const hal_index_t parentIndex = topSegment->getParentIndex();
  if (parentGenome != NULL && parentIndex != NULL_INDEX)
  {
    if (parentIndex >= (hal_index_t)parentGenome->getNumBottomSegments())
    {
        throw hal_exception("Parent index " + std::to_string(parentIndex) + " of segment "
                            + std::to_string(topSegment->getArrayIndex()) + " out of range in genome "
                            + parentGenome->getName());
    }
    BottomSegmentIteratorPtr bottomSegmentIterator = 
       parentGenome->getBottomSegmentIterator(parentIndex);
    const BottomSegment* parentSegment = 
       bottomSegmentIterator->getBottomSegment();
    if (topSegment->getLength() != parentSegment->getLength())
    {
        throw hal_exception("Parent length of segment " + std::to_string(topSegment->getArrayIndex()) 
                            + " in genome " + genome->getName() + " has length "
                            + std::to_string(parentSegment->getLength()) + " which does not match "
                            + std::to_string(topSegment->getLength()));
    }
  }

  const hal_index_t parseIndex = topSegment->getBottomParseIndex();
  if (parseIndex == NULL_INDEX)
  {
    if (genome->getNumChildren() != 0)
    {
        throw hal_exception("Top Segment " + std::to_string(topSegment->getArrayIndex()) + " in genome "
                            + genome->getName() + " has null parse index");
    }
  }
  else
  {
    if (parseIndex >= (hal_index_t)genome->getNumBottomSegments())
    {
        throw hal_exception("Top Segment " + std::to_string(topSegment->getArrayIndex()) + " in genome "
                            + genome->getName() + " has parse index " + std::to_string(parseIndex)
                            + " which is out of range since genome has " 
                            + std::to_string(genome->getNumBottomSegments()) + " bottom segments");
    }
    hal_offset_t parseOffset = topSegment->getBottomParseOffset();
    BottomSegmentIteratorPtr bottomSegmentIterator =
       genome->getBottomSegmentIterator(parseIndex);
    const BottomSegment* parseSegment = 
       bottomSegmentIterator->getBottomSegment();
    if (parseOffset >= parseSegment->getLength())
    {
        throw hal_exception("Top Segment " + std::to_string(topSegment->getArrayIndex()) + " in genome "
                          + genome->getName() + " has parse offset out of range");
    }
    if ((hal_index_t)parseOffset + parseSegment->getStartPosition() != 
        topSegment->getStartPosition())
    {
      throw hal_exception("parse index broken in top segment in genome " +
                          genome->getName());
                          
    }
  }

  const hal_index_t paralogyIndex = topSegment->getNextParalogyIndex();
  if (paralogyIndex != NULL_INDEX)
  {
    TopSegmentIteratorPtr pti = 
       genome->getTopSegmentIterator(paralogyIndex);
    if (pti->getTopSegment()->getParentIndex() != topSegment->getParentIndex())
    {
        throw hal_exception("Top segment " + std::to_string(topSegment->getArrayIndex())
                            + " has parent index " + std::to_string(topSegment->getParentIndex()) + ", but next paraglog " 
                            + std::to_string(topSegment->getNextParalogyIndex()) + " has parent Index " 
                            + std::to_string(pti->getTopSegment()->getParentIndex()) 
                            + ". Paralogous top segments must share same parent.");
    }
    if (paralogyIndex == topSegment->getArrayIndex())
    {
        throw hal_exception("Top segment " + std::to_string(topSegment->getArrayIndex())
                            + " has paralogy index " + std::to_string(topSegment->getNextParalogyIndex())
                            + " which isn't allowed");
    }
  }
}

void hal::validateSequence(const Sequence* sequence)
{
  // Verify that the DNA sequence doesn't contain funny characters
  DnaIteratorPtr dnaIt = sequence->getDnaIterator();
  hal_size_t length = sequence->getSequenceLength();
  if (sequence->getGenome()->containsDNAArray() == true)
  {
    for (hal_size_t i = 0; i < length; ++i)
    {
      char c = dnaIt->getBase();
      if (isNucleotide(c) == false)
      {
        throw hal_exception( "Non-nucleotide character discoverd at position " 
                             + std::to_string(i) + " of sequence " + sequence->getName() + ": " + c);
      }
    }
  }

  // Check the top segments
  if (sequence->getGenome()->getParent() != NULL)
  {
    hal_size_t totalTopLength = 0;
    TopSegmentIteratorPtr topIt = sequence->getTopSegmentIterator();
    hal_size_t numTopSegments = sequence->getNumTopSegments();
    for (hal_size_t i = 0; i < numTopSegments; ++i)
    {
      const TopSegment* topSegment = topIt->getTopSegment();
      validateTopSegment(topSegment);
      totalTopLength += topSegment->getLength();
      topIt->toRight();
    }
    if (totalTopLength != length)
    {
        throw hal_exception("Sequence " + sequence->getName() + " has length " + std::to_string(length )
                            + " but its top segments add up to " + std::to_string(totalTopLength));
    }
  }

  // Check the bottom segments
  if (sequence->getGenome()->getNumChildren() > 0)
  {
    hal_size_t totalBottomLength = 0;
    BottomSegmentIteratorPtr bottomIt = sequence->getBottomSegmentIterator();
    hal_size_t numBottomSegments = sequence->getNumBottomSegments();
    for (hal_size_t i = 0; i < numBottomSegments; ++i)
    {
      const BottomSegment* bottomSegment = bottomIt->getBottomSegment();
      validateBottomSegment(bottomSegment);
      totalBottomLength += bottomSegment->getLength();
      bottomIt->toRight();
    }
    if (totalBottomLength != length)
    {
        throw hal_exception("Sequence " + sequence->getName() + " has length " + std::to_string(length)
                            + " but its bottom segments add up to " + std::to_string(totalBottomLength));
    }
  }
}

void hal::validateDuplications(const Genome* genome)
{
  const Genome* parent = genome->getParent();
  if (parent == NULL)
  {
    return;
  }
  TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
  vector<unsigned char> pcount(parent->getNumBottomSegments(), 0);
  for (; not topIt->atEnd(); topIt->toRight())
  {
    if (topIt->tseg()->hasParent())
    {
      if (pcount[topIt->getTopSegment()->getParentIndex()] < 250)
      {
        ++pcount[topIt->getTopSegment()->getParentIndex()];
      }
    }
  }
  for (topIt = genome->getTopSegmentIterator(); not topIt->atEnd(); topIt->toRight())
  {
    if (topIt->tseg()->hasParent())
    {
      size_t count = pcount[topIt->getTopSegment()->getParentIndex()];
      assert(count > 0);
      {
        if (topIt->tseg()->hasNextParalogy() == false && count > 1)
        {
            throw hal_exception("Top Segment " + std::to_string(topIt->getTopSegment()->getArrayIndex())
                                + " in genome " + genome->getName() + " is not marked as a"
                                + " duplication but it shares its parent " 
                                + std::to_string(topIt->getTopSegment()->getArrayIndex()) + " with at least " 
                                                 + std::to_string(count - 1) + " other segments in the same genome");
        }  
      }
    }
  }
}

void hal::validateGenome(const Genome* genome)
{
  // first we check the sequence coverage
  hal_size_t totalTop = 0;
  hal_size_t totalBottom = 0;
  hal_size_t totalLength = 0;
  
  SequenceIteratorPtr seqIt = genome->getSequenceIterator();
  for (; not seqIt->atEnd(); seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    validateSequence(sequence);

    totalTop += sequence->getNumTopSegments();
    totalBottom += sequence->getNumBottomSegments();
    totalLength += sequence->getSequenceLength();

    // make sure it doesn't overlap any other sequences;
    if (sequence->getSequenceLength() > 0)
    {
      const Sequence* s1 =
         genome->getSequenceBySite(sequence->getStartPosition());

      if (s1 == NULL || s1->getName() != sequence->getName())
      {
          throw hal_exception("Sequence " + sequence->getName() + " has a bad overlap in "
                              + genome->getName());
      }
      const Sequence* s2 = 
         genome->getSequenceBySite(sequence->getStartPosition() +
                                   sequence->getSequenceLength() - 1);
      if (s2 == NULL || s2->getName() != sequence->getName())
      {
        throw hal_exception("Sequence " + sequence->getName() + " has a bad overlap in "
                             + genome->getName());
      }
    }
  }

  hal_size_t genomeLength = genome->getSequenceLength();
  hal_size_t genomeTop = genome->getNumTopSegments();
  hal_size_t genomeBottom = genome->getNumBottomSegments();

  if (genomeLength != totalLength)
  {
      throw hal_exception("Problem: genome has length " + std::to_string(genomeLength)
                          + ", however sequences total " + std::to_string(totalLength));
  }
  if (genomeTop != totalTop)
  {
      throw hal_exception("Problem: genome has " + std::to_string(genomeTop) + " top segments but "
                          + "sequences have " + std::to_string(totalTop) + " top segments");
  }
  if (genomeBottom != totalBottom)
  {
      throw hal_exception("Problem: genome has " + std::to_string(genomeBottom) + " bottom segments but "
                          + "sequences have " + std::to_string(totalBottom) + " bottom segments");
  }

  if (genomeLength > 0 && genomeTop == 0 && genomeBottom == 0)
  {
    throw hal_exception("Problem: genome " + genome->getName() + " has length " 
                        + std::to_string(genomeLength) + "but no segments");
  }
  
  validateDuplications(genome);
}

void hal::validateAlignment(const Alignment* alignment)
{
  deque<string> bfQueue;
  bfQueue.push_back(alignment->getRootName());
  while (bfQueue.empty() == false)
  {
    string name = bfQueue.back();
    bfQueue.pop_back();
    if (name.empty() == false)
    {
        const Genome* genome = alignment->openGenome(name);
      if (genome == NULL)
      {
        throw hal_exception("Failure to open genome " + name);
      }
      validateGenome(genome);
      vector<string> childNames = alignment->getChildNames(name);
      for (size_t i = 0; i < childNames.size(); ++i)
      {
        bfQueue.push_front(childNames[i]);
      }
    }
  }
}
