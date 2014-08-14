/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <sstream>
#include <deque>
#include <vector>
#include <iostream>
#include "halValidate.h"
#include "hal.h"

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
    stringstream ss;
    ss << "Bottom segment out of range " << index << " in genome "
       << genome->getName();
    throw hal_exception(ss.str());
  }
  
  if (bottomSegment->getLength() < 1)
  {
    stringstream ss;
    ss << "Bottom segment " << index  << " in genome " << genome->getName()
       << " has length 0 which is not currently supported";
    throw hal_exception(ss.str());
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
        stringstream ss;
        ss << "Child " << child << " index " <<childIndex << " of segment "
           << bottomSegment->getArrayIndex() << " out of range in genome "
           << childGenome->getName();
        throw hal_exception(ss.str());
      }
      TopSegmentIteratorConstPtr topSegmentIteratr = 
         childGenome->getTopSegmentIterator(childIndex);
      const TopSegment* childSegment = topSegmentIteratr->getTopSegment();
      if (childSegment->getLength() != bottomSegment->getLength())
      {
        stringstream ss;
        ss << "Child " << child << " with index " 
           << childSegment->getArrayIndex()
           << " and start position " << childSegment->getStartPosition() 
           << " and sequence " << childSegment->getSequence()->getName()
           << " has length " << childSegment->getLength()
           << " but parent with index " << bottomSegment->getArrayIndex() 
           << " and start position " << bottomSegment->getStartPosition()
           << " in sequence " << bottomSegment->getSequence()->getName() 
           << " has length " << bottomSegment->getLength();
        throw hal_exception(ss.str());
      }
      if (childSegment->getNextParalogyIndex() == NULL_INDEX &&
          childSegment->getParentIndex() != bottomSegment->getArrayIndex())
      {
        stringstream ss;
        ss << "Parent / child index mismatch:\n" 
           << genome->getName() << "[" << bottomSegment->getArrayIndex() << "]"
           << " links to " << childGenome->getName() << "[" << childIndex 
           << "] but \n"
           << childGenome->getName() << "[" << childSegment->getArrayIndex() 
           << "] links to " << genome->getName() << "[" 
           << childSegment->getParentIndex() << "]";
        throw hal_exception(ss.str());
      }
      if (childSegment->getParentReversed() != 
          bottomSegment->getChildReversed(child))
      {
        stringstream ss;
        ss << "parent / child reversal mismatch (parent=" <<
          genome->getName() << " parentSegNum=" << bottomSegment->getArrayIndex() << " child=" <<
          childGenome->getName() << " childSegNum=" << childSegment->getArrayIndex() << ")";
          throw hal_exception(ss.str());
      }
    }
  }

  const hal_index_t parseIndex = bottomSegment->getTopParseIndex();
  if (parseIndex == NULL_INDEX)
  {
    if (genome->getParent() != NULL)
    {
      stringstream ss;
      ss << "Bottom segment " << bottomSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has null parse index";
      throw hal_exception(ss.str());
    }
  }
  else
  {
    if (parseIndex >= (hal_index_t)genome->getNumTopSegments())
    {
      stringstream ss;
      ss << "BottomSegment " << bottomSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has parse index " << parseIndex 
         << " greater than the number of top segments, " 
         << (hal_index_t)genome->getNumTopSegments();
      throw hal_exception(ss.str());
    }
    TopSegmentIteratorConstPtr parseIterator = 
       genome->getTopSegmentIterator(parseIndex);
    const TopSegment* parseSegment = parseIterator->getTopSegment();
    hal_offset_t parseOffset = bottomSegment->getTopParseOffset();
    if (parseOffset >= parseSegment->getLength())
    {
      stringstream ss;
      ss << "BottomSegment " << bottomSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has parse offset, " << parseOffset 
         << ", greater than the length of the segment, " 
         << parseSegment->getLength();
      throw hal_exception(ss.str());
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
    stringstream ss;
    ss << "Segment out of range " << index << " in genome "
       << genome->getName();
    throw hal_exception(ss.str());
  }

  if (topSegment->getLength() < 1)
  {
    stringstream ss;
    ss << "Top segment " << index  << " in genome " << genome->getName()
       << " has length 0 which is not currently supported";
    throw hal_exception(ss.str());
  }

  const Genome* parentGenome = genome->getParent();
  const hal_index_t parentIndex = topSegment->getParentIndex();
  if (parentGenome != NULL && parentIndex != NULL_INDEX)
  {
    if (parentIndex >= (hal_index_t)parentGenome->getNumBottomSegments())
    {
      stringstream ss;
      ss << "Parent index " << parentIndex << " of segment "
         << topSegment->getArrayIndex() << " out of range in genome "
         << parentGenome->getName();
      throw hal_exception(ss.str());
    }
    BottomSegmentIteratorConstPtr bottomSegmentIterator = 
       parentGenome->getBottomSegmentIterator(parentIndex);
    const BottomSegment* parentSegment = 
       bottomSegmentIterator->getBottomSegment();
    if (topSegment->getLength() != parentSegment->getLength())
    {
      stringstream ss;
      ss << "Parent length of segment " << topSegment->getArrayIndex() 
         << " in genome " << genome->getName() << " has length "
         << parentSegment->getLength() << " which does not match "
         << topSegment->getLength();
      throw hal_exception(ss.str());
    }
  }

  const hal_index_t parseIndex = topSegment->getBottomParseIndex();
  if (parseIndex == NULL_INDEX)
  {
    if (genome->getNumChildren() != 0)
    {
      stringstream ss;
      ss << "Top Segment " << topSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has null parse index";
      throw hal_exception(ss.str());
    }
  }
  else
  {
    if (parseIndex >= (hal_index_t)genome->getNumBottomSegments())
    {
      stringstream ss;
      ss << "Top Segment " << topSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has parse index " << parseIndex 
         << " which is out of range since genome has " 
         << genome->getNumBottomSegments() << " bottom segments";
      throw hal_exception(ss.str());
    }
    hal_offset_t parseOffset = topSegment->getBottomParseOffset();
    BottomSegmentIteratorConstPtr bottomSegmentIterator =
       genome->getBottomSegmentIterator(parseIndex);
    const BottomSegment* parseSegment = 
       bottomSegmentIterator->getBottomSegment();
    if (parseOffset >= parseSegment->getLength())
    {
      stringstream ss;
      ss << "Top Segment " << topSegment->getArrayIndex() << " in genome "
         << genome->getName() << " has parse offset out of range";
      throw hal_exception(ss.str());
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
    TopSegmentIteratorConstPtr pti = 
       genome->getTopSegmentIterator(paralogyIndex);
    if (pti->getTopSegment()->getParentIndex() != topSegment->getParentIndex())
    {
      stringstream ss;
      ss << "Top segment " << topSegment->getArrayIndex() 
         << " has parent index "
         << topSegment->getParentIndex() << ", but next paraglog " 
         << topSegment->getNextParalogyIndex() << " has parent Index " 
         << pti->getTopSegment()->getParentIndex() 
         << ". Paralogous top segments must share same parent.";
      throw hal_exception(ss.str());
    }
    if (paralogyIndex == topSegment->getArrayIndex())
    {
      stringstream ss;
      ss << "Top segment " << topSegment->getArrayIndex() 
         << " has paralogy index " << topSegment->getNextParalogyIndex()
         << " which isn't allowed";
      throw hal_exception(ss.str());
    }
  }
}

void hal::validateSequence(const Sequence* sequence)
{
  // Verify that the DNA sequence doesn't contain funny characters
  DNAIteratorConstPtr dnaIt = sequence->getDNAIterator();
  hal_size_t length = sequence->getSequenceLength();
  if (sequence->getGenome()->containsDNAArray() == true)
  {
    for (hal_size_t i = 0; i < length; ++i)
    {
      char c = dnaIt->getChar();
      if (isNucleotide(c) == false)
      {
        stringstream ss;
        ss << "Non-nucleotide character discoverd at position " 
           << i << " of sequence " << sequence->getName() << ": " << c;
        throw hal_exception(ss.str());
      }
    }
  }

  // Check the top segments
  if (sequence->getGenome()->getParent() != NULL)
  {
    hal_size_t totalTopLength = 0;
    TopSegmentIteratorConstPtr topIt = sequence->getTopSegmentIterator();
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
      stringstream ss;
      ss << "Sequence " << sequence->getName() << " has length " << length 
         << " but its top segments add up to " << totalTopLength;
      throw hal_exception(ss.str());
    }
  }

  // Check the bottom segments
  if (sequence->getGenome()->getNumChildren() > 0)
  {
    hal_size_t totalBottomLength = 0;
    BottomSegmentIteratorConstPtr bottomIt = 
       sequence->getBottomSegmentIterator();
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
      stringstream ss;
      ss << "Sequence " << sequence->getName() << " has length " << length 
         << " but its bottom segments add up to " << totalBottomLength;
      throw hal_exception(ss.str());
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
  TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
  TopSegmentIteratorConstPtr endIt = genome->getTopSegmentEndIterator();
  vector<unsigned char> pcount(parent->getNumBottomSegments(), 0);
  for (; topIt != endIt; topIt->toRight())
  {
    if (topIt->hasParent())
    {
      if (pcount[topIt->getTopSegment()->getParentIndex()] < 250)
      {
        ++pcount[topIt->getTopSegment()->getParentIndex()];
      }
    }
  }
  for (topIt = genome->getTopSegmentIterator(); topIt != endIt; topIt->toRight())
  {
    if (topIt->hasParent())
    {
      size_t count = pcount[topIt->getTopSegment()->getParentIndex()];
      assert(count > 0);
      {
        if (topIt->hasNextParalogy() == false && count > 1)
        {
          stringstream ss;
          ss << "Top Segment " << topIt->getTopSegment()->getArrayIndex()
             << " in genome " << genome->getName() << " is not marked as a"
             << " duplication but it shares its parent " 
             << topIt->getTopSegment()->getArrayIndex() << " with at least " 
             << count - 1 << " other segments in the same genome";
          throw hal_exception(ss.str());
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
  
  SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
  SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
  for (; seqIt != seqEnd; seqIt->toNext())
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
        stringstream ss;
        ss << "Sequence " << sequence->getName() << " has a bad overlap in "
           << genome->getName();
        throw hal_exception(ss.str());
      }
      const Sequence* s2 = 
         genome->getSequenceBySite(sequence->getStartPosition() +
                                   sequence->getSequenceLength() - 1);
      if (s2 == NULL || s2->getName() != sequence->getName())
      {
        stringstream ss;
        ss << "Sequence " << sequence->getName() << " has a bad overlap in "
           << genome->getName();
        throw hal_exception(ss.str());
      }
    }
  }

  hal_size_t genomeLength = genome->getSequenceLength();
  hal_size_t genomeTop = genome->getNumTopSegments();
  hal_size_t genomeBottom = genome->getNumBottomSegments();

  if (genomeLength != totalLength)
  {
    stringstream ss;
    ss << "Problem: genome has length " << genomeLength 
       << "But sequences total " << totalLength;
    throw hal_exception(ss.str());
  }
  if (genomeTop != totalTop)
  {
    stringstream ss;
    ss << "Problem: genome has " << genomeTop << " top segments but "
       << "sequences have " << totalTop << " top segments";
    throw ss.str();
  }
  if (genomeBottom != totalBottom)
  {
    stringstream ss;
    ss << "Problem: genome has " << genomeBottom << " bottom segments but "
       << "sequences have " << totalBottom << " bottom segments";
    throw hal_exception(ss.str());
  }

  if (genomeLength > 0 && genomeTop == 0 && genomeBottom == 0)
  {
    stringstream ss;
    ss << "Problem: genome " << genome->getName() << " has length " 
       << genomeLength << "but no segments";
    throw hal_exception(ss.str());
  }
  
  validateDuplications(genome);
}

void hal::validateAlignment(AlignmentConstPtr alignment)
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
