#include <sstream>
#include <assert.h>
#include <map>
#include <iostream>
#include "halGenome.h"
#include "halAlignment.h"
#include "halBottomSegmentIterator.h"
#include "halDNAIterator.h"
#include "halMetaData.h"
#include "halSegmentIterator.h"
#include "halSequenceIterator.h"
#include "halTopSegmentIterator.h"
using namespace std;

namespace hal {
void Genome::copy(Genome *dest) const
{
  copyDimensions(dest);
  copySequence(dest);
  copyTopSegments(dest);
  copyBottomSegments(dest);
  copyMetadata(dest);
}

void Genome::copyDimensions(Genome *dest) const
{
  vector<Sequence::Info> dimensions;
  const Alignment *inAlignment = getAlignment();
  SequenceIteratorConstPtr seqIt = getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = getSequenceEndIterator();

  bool root = inAlignment->getParentName(getName()).empty();
  bool leaf = inAlignment->getChildNames(getName()).empty();     
  
  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    Sequence::Info info(sequence->getName(),
                        sequence->getSequenceLength(),
                        root ? 0 : sequence->getNumTopSegments(),
                        leaf ? 0 : sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
  dest->setDimensions(dimensions);
}

void Genome::copyTopDimensions(Genome *dest) const
{
  vector<Sequence::UpdateInfo> dimensions;
  SequenceIteratorConstPtr seqIt = getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = getSequenceEndIterator();
  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    if (sequence->getSequenceLength() == 0 &&
        dest->getSequence(sequence->getName()) == NULL) {
      // progressiveCactus creates 0-length sequences in ancestors,
      // which are not usually extractable and aren't important
      continue;
    }
    Sequence::UpdateInfo info(sequence->getName(),
                              sequence->getNumTopSegments());
    dimensions.push_back(info);
  }
  dest->updateTopDimensions(dimensions);
}

void Genome::copyBottomDimensions(Genome *dest) const
{
  vector<Sequence::UpdateInfo> dimensions;
  SequenceIteratorConstPtr seqIt = getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = getSequenceEndIterator();

  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    if (sequence->getSequenceLength() == 0 &&
        dest->getSequence(sequence->getName()) == NULL) {
      // progressiveCactus creates 0-length sequences in ancestors,
      // which are not usually extractable and aren't important
      continue;
    }
    Sequence::UpdateInfo info(sequence->getName(),
                              sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
  dest->updateBottomDimensions(dimensions);
}

void Genome::copyTopSegments(Genome *dest) const
{
  const Genome *inParent = getParent();
  const Genome *outParent = dest->getParent();

  TopSegmentIteratorConstPtr inTop = getTopSegmentIterator();
  TopSegmentIteratorPtr outTop = dest->getTopSegmentIterator();
  hal_size_t n = dest->getNumTopSegments();
  assert(n == 0 || n == getNumTopSegments());

  if (n == 0) {
    // Nothing to do if there are no top segments.
    return;
  }

  BottomSegmentIteratorConstPtr inParentBottomSegIt = inParent->getBottomSegmentIterator();
  BottomSegmentIteratorConstPtr outParentBottomSegIt = outParent->getBottomSegmentIterator();


  // Go through each sequence in this genome, find the matching
  // sequence in the dest genome, then copy over the segments for each
  // sequence.
  SequenceIteratorConstPtr seqIt = getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = getSequenceEndIterator();

  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence *inSeq = seqIt->getSequence();
    const Sequence *outSeq = dest->getSequence(inSeq->getName());
    TopSegmentIteratorPtr inTop = inSeq->getTopSegmentIterator();
    TopSegmentIteratorPtr outTop = outSeq->getTopSegmentIterator();

    if (inSeq->getName() != outSeq->getName()) {
      stringstream ss;
      ss << "When copying top segments: segment #" << inTop->getArrayIndex() << " of source genome is from sequence " << inTop->getSequence()->getName() << ", but segment #" << outTop->getArrayIndex() << " is from sequence " << outTop->getSequence()->getName();
      throw hal_exception(ss.str());
    }

    if (inSeq->getNumTopSegments() != outSeq->getNumTopSegments()) {
      stringstream ss;
      ss << "When copying top segments: sequence " << inSeq->getName() << " has " << inSeq->getNumTopSegments() << " in genome " << getName() << ", while it has " << outSeq->getNumTopSegments() << " in genome " << dest->getName();
      throw hal_exception(ss.str());
    }

    hal_index_t inSegmentEnd = inSeq->getTopSegmentArrayIndex() + inSeq->getNumTopSegments();
    for (; (hal_size_t)inTop->getArrayIndex() < inSegmentEnd; inTop->toRight(),
             outTop->toRight())
    {
        hal_index_t outStartPosition = inTop->getStartPosition() - inSeq->getStartPosition() + outSeq->getStartPosition();
        outTop->setCoordinates(outStartPosition, inTop->getLength());
        outTop->setParentIndex(inTop->getParentIndex());
        outTop->setParentReversed(inTop->getParentReversed());
        outTop->setBottomParseIndex(inTop->getBottomParseIndex());
        // Figure out next paralogy index in new sequence order.
        hal_index_t inParalogyIndex = inTop->getNextParalogyIndex();
        hal_index_t outParalogyIndex;
        if (inParalogyIndex != NULL_INDEX) {
            TopSegmentIteratorPtr inParalogyIt = getTopSegmentIterator(inParalogyIndex);
            Sequence *inParalogySeq = inParalogyIt->getSequence();
            Sequence *outParalogySeq = dest->getSequence(inParalogySeq->getName());
            outParalogyIndex = inParalogyIndex - inParalogySeq->getTopSegmentArrayIndex() + outParalogySeq->getTopSegmentArrayIndex();
        } else {
            outParalogyIndex = NULL_INDEX;
        }
        outTop->setNextParalogyIndex(outParalogyIndex);

        // Check that the sequences from the bottom segments we point to are the same. If not, correct the indices so that they are.
        if (inTop->getParentIndex() != NULL_INDEX) {
            inParentBottomSegIt->toParent(inTop);

            const Sequence *inParentSequence = inParentBottomSegIt->getSequence();

            const Sequence *outParentSequence = outParent->getSequence(inParentSequence->getName());

            hal_index_t inParentSegmentOffset = inTop->getParentIndex() - inParentSequence->getBottomSegmentArrayIndex();
            hal_index_t outParentSegmentIndex = inParentSegmentOffset + outParentSequence->getBottomSegmentArrayIndex();

            outTop->setParentIndex(outParentSegmentIndex);
        } else {
            outTop->setParentIndex(NULL_INDEX);
        }
    }
  }
}

void Genome::copyBottomSegments(Genome *dest) const
{
  assert(getNumBottomSegments() == dest->getNumBottomSegments());
  hal_size_t inNc = getNumChildren();
  hal_size_t outNc = dest->getNumChildren();
  // The child indices aren't consistent across files--make sure each bottom
  // segment points to the correct children
  vector<string> inChildNames;
  vector<string> outChildNames;
  for (hal_size_t inChild = 0; inChild < inNc; ++inChild)
  {
    inChildNames.push_back(getChild(inChild)->getName());
  }
  for (hal_size_t outChild = 0; outChild < outNc; ++outChild)
  {
    outChildNames.push_back(dest->getChild(outChild)->getName());
  }
  map<hal_size_t, hal_size_t> inChildToOutChild;
  for (hal_size_t inChild = 0; inChild < inNc; inChild++)
  {
    hal_size_t outChild;
    for (outChild = 0; outChild < outNc; outChild++)
    {
      if (inChildNames[inChild] == outChildNames[outChild])
      {
        inChildToOutChild[inChild] = outChild;
        break;
      }
    }
    if (outChild == outNc)
    {
      inChildToOutChild[inChild] = outNc;
    }
  }

  // Go through each sequence in this genome, find the matching
  // sequence in the dest genome, then copy over the segments for each
  // sequence.
  SequenceIteratorConstPtr seqIt = getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = getSequenceEndIterator();

  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence *inSeq = seqIt->getSequence();
    const Sequence *outSeq = dest->getSequence(inSeq->getName());
    BottomSegmentIteratorPtr inBot = inSeq->getBottomSegmentIterator();
    BottomSegmentIteratorPtr outBot = outSeq->getBottomSegmentIterator();

    if (inSeq->getName() != outSeq->getName()) {
      // This check is important enough that it can't be an assert.
      stringstream ss;
      ss << "When copying bottom segments: segment #" << inBot->getArrayIndex() << " of source genome is from sequence " << inBot->getSequence()->getName() << ", but segment #" << outBot->getArrayIndex() << " is from sequence " << outBot->getSequence()->getName();
      throw hal_exception(ss.str());
    }

    if (inSeq->getNumBottomSegments() != outSeq->getNumBottomSegments()) {
      stringstream ss;
      ss << "When copying bottom segments: sequence " << inSeq->getName() << " has " << inSeq->getNumBottomSegments() << " in genome " << getName() << ", while it has " << outSeq->getNumBottomSegments() << " in genome " << dest->getName();
      throw hal_exception(ss.str());      
    }

    hal_index_t inSegmentEnd = inSeq->getBottomSegmentArrayIndex() + inSeq->getNumBottomSegments();
    for (; inBot->getArrayIndex() < inSegmentEnd; inBot->toRight(),
           outBot->toRight())
    {
      hal_index_t outStartPosition = inBot->getStartPosition() - inSeq->getStartPosition() + outSeq->getStartPosition();

      if (dest->getSequenceBySite(outStartPosition) != outSeq) {
        stringstream ss;
        ss << "When copying bottom segments from " << getName() << " to " << dest->getName() << ": expected destination sequence " << outSeq->getName() << " for segment # " << inBot->getArrayIndex() << " but got " << dest->getSequenceBySite(outStartPosition)->getName();
        throw hal_exception(ss.str());
      }
      outBot->setCoordinates(outStartPosition, inBot->getLength());
      for(hal_size_t inChild = 0; inChild < inNc; inChild++) {
        hal_size_t outChild = inChildToOutChild[inChild];
        if (outChild != outNc) {
            // Adjust top segment we point to so that it points to the
            // same sequence, even if the sequence ordering has changed.
            hal_index_t childIndex = inBot->getChildIndex(inChild);
            Genome *outChildGenome = dest->getChild(outChild);
            if (childIndex != NULL_INDEX) {
                TopSegmentIteratorConstPtr inIt = getTopSegmentIterator();
                inIt->toChild(inBot, inChild);

                const Sequence *inChildSequence = inIt->getSequence();
                const Sequence *outChildSequence = outChildGenome->getSequence(inChildSequence->getName());

                hal_index_t inChildSegmentOffset = inIt->getArrayIndex() - inChildSequence->getTopSegmentArrayIndex();
                childIndex = inChildSegmentOffset + outChildSequence->getTopSegmentArrayIndex();
            }

            outBot->setChildIndex(outChild, childIndex);
            outBot->setChildReversed(outChild, inBot->getChildReversed(inChild));
        }
      }
      outBot->setTopParseIndex(inBot->getTopParseIndex());
    }
  }
}


void Genome::copySequence(Genome *dest) const
{
  DNAIteratorConstPtr inDna = getDNAIterator();
  DNAIteratorPtr outDna = dest->getDNAIterator();
  hal_size_t n = getSequenceLength();
  assert(n == dest->getSequenceLength());
  for (; (hal_size_t)inDna->getArrayIndex() < n; inDna->toRight(), 
         outDna->toRight())
  {
    outDna->setChar(inDna->getChar());
  }
}

void Genome::copyMetadata(Genome *dest) const
{
  const map<string, string>& meta = getMetaData()->getMap();
  map<string, string>::const_iterator i = meta.begin();
  for (; i != meta.end(); ++i)
  {
    dest->getMetaData()->set(i->first, i->second);
  }
}

void Genome::fixParseInfo()
{
  if (getParent() == NULL || getNumChildren() == 0)
  {
    return;
  }
  
  // copied from CactusHalConverter::updateRootParseInfo() in
  // cactus2hal/src/cactusHalConverter.cpp 
  BottomSegmentIteratorPtr bottomIterator = 
    getBottomSegmentIterator();
  TopSegmentIteratorPtr topIterator = getTopSegmentIterator();
  BottomSegmentIteratorConstPtr bend = getBottomSegmentEndIterator();
  TopSegmentIteratorConstPtr tend = getTopSegmentEndIterator();
  int top = 0, bot = 0;
  while (bottomIterator != bend && topIterator != tend)
  {
    bool bright = false;
    bool tright = false;
    BottomSegment* bseg = bottomIterator->getBottomSegment();
    TopSegment* tseg = topIterator->getTopSegment();
    hal_index_t bstart = bseg->getStartPosition();
    hal_index_t bendidx = bstart + (hal_index_t)bseg->getLength();
    hal_index_t tstart = tseg->getStartPosition();
    hal_index_t tendidx = tstart + (hal_index_t)tseg->getLength();

    if (bstart >= tstart && bstart < tendidx)
    {
      bseg->setTopParseIndex(tseg->getArrayIndex());
    }
    if (bendidx <= tendidx || bstart == bendidx)
    {
      bright = true;
    }
        
    if (tstart >= bstart && tstart < bendidx)
    {
      tseg->setBottomParseIndex(bseg->getArrayIndex());
    }
    if (tendidx <= bendidx || tstart == tendidx)
    {
      tright = true;
    }

    assert(bright || tright);
    if (bright == true)
    {
      bot += 1;
      bottomIterator->toRight();
    }
    if (tright == true)
    {
      top += 1;
      topIterator->toRight();
    }
  }
}
} // namespace hal
