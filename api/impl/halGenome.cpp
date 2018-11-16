#include <cassert>
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
  SequenceIteratorPtr seqIt = getSequenceIterator();
  SequenceIteratorPtr seqEndIt = getSequenceEndIterator();

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
  SequenceIteratorPtr seqIt = getSequenceIterator();
  SequenceIteratorPtr seqEndIt = getSequenceEndIterator();
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
  SequenceIteratorPtr seqIt = getSequenceIterator();
  SequenceIteratorPtr seqEndIt = getSequenceEndIterator();

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

  TopSegmentIteratorPtr inTop = getTopSegmentIterator();
  TopSegmentIteratorPtr outTop = dest->getTopSegmentIterator();
  hal_size_t n = dest->getNumTopSegments();
  assert(n == 0 || n == getNumTopSegments());

  if (n == 0) {
    // Nothing to do if there are no top segments.
    return;
  }

  BottomSegmentIteratorPtr inParentBottomSegIt = inParent->getBottomSegmentIterator();
  BottomSegmentIteratorPtr outParentBottomSegIt = outParent->getBottomSegmentIterator();

  for (; (hal_size_t)inTop->getArrayIndex() < n; inTop->toRight(),
         outTop->toRight())
  {
    hal_index_t genomePos = inTop->getStartPosition();
    assert(genomePos != NULL_INDEX);
    string inSeqName = getSequenceBySite(genomePos)->getName();
    string outSeqName = dest->getSequenceBySite(genomePos)->getName();
#if NDEBUG
    if (inSeqName != outSeqName) {
        throw hal_exception("When copying top segments from " + getName() + " to " + dest->getName() + ": sequence " + inSeqName + " != " + outSeqName + " at site " + std::to_string(genomePos));
    }
#endif

    outTop->setCoordinates(inTop->getStartPosition(), inTop->getLength());
    outTop->setParentIndex(inTop->getParentIndex());
    outTop->setParentReversed(inTop->getParentReversed());
    outTop->setBottomParseIndex(inTop->getBottomParseIndex());
    outTop->setNextParalogyIndex(inTop->getNextParalogyIndex());

    // Check that the sequences from the bottom segments we point to are the same. If not, correct the indices so that they are.
    if (inTop->getParentIndex() != NULL_INDEX) {
      inParentBottomSegIt->toParent(inTop);

      const Sequence *inParentSequence = inParentBottomSegIt->getSequence();

      const Sequence *outParentSequence = outParent->getSequence(inParentSequence->getName());

      hal_index_t inParentSegmentOffset = inTop->getParentIndex() - inParentSequence->getBottomSegmentArrayIndex();
      hal_index_t outParentSegmentIndex = inParentSegmentOffset + outParentSequence->getBottomSegmentArrayIndex();

      outTop->setParentIndex(outParentSegmentIndex);
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
  SequenceIteratorPtr seqIt = getSequenceIterator();
  SequenceIteratorPtr seqEndIt = getSequenceEndIterator();

  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence *inSeq = seqIt->getSequence();
    const Sequence *outSeq = dest->getSequence(inSeq->getName());
    BottomSegmentIteratorPtr inBotSegIt = const_pointer_cast<BottomSegmentIterator>(inSeq->getBottomSegmentIterator());
    BottomSegmentIteratorPtr outBotSegIt = const_pointer_cast<BottomSegmentIterator>(outSeq->getBottomSegmentIterator());

    if (inSeq->getName() != outSeq->getName()) {
      // This check is important enough that it can't be an assert.
      throw hal_exception("When copying bottom segments: segment #" + std::to_string(inBotSegIt->getArrayIndex()) + " of source genome is from sequence " + inBotSegIt->getSequence()->getName() + ", but segment #" + std::to_string(outBotSegIt->getArrayIndex()) + " is from sequence " + outBotSegIt->getSequence()->getName());
    }

    if (inSeq->getNumBottomSegments() != outSeq->getNumBottomSegments()) {
        throw hal_exception("When copying bottom segments: sequence " + inSeq->getName() + " has " + std::to_string(inSeq->getNumBottomSegments()) + " in genome " + getName() + ", while it has " + std::to_string(outSeq->getNumBottomSegments()) + " in genome " + dest->getName());
    }

    hal_index_t inSegmentEnd = inSeq->getBottomSegmentArrayIndex() + inSeq->getNumBottomSegments();
    for (; inBotSegIt->getArrayIndex() < inSegmentEnd; inBotSegIt->toRight(),
           outBotSegIt->toRight())
    {
      hal_index_t outStartPosition = inBotSegIt->getStartPosition() - inSeq->getStartPosition() + outSeq->getStartPosition();


      if (dest->getSequenceBySite(outStartPosition) != outSeq) {
          throw hal_exception("When copying bottom segments from " + getName() + " to " + dest->getName() + ": expected destination sequence " + outSeq->getName() + " for segment # " + std::to_string(inBotSegIt->getArrayIndex()) + " but got " + dest->getSequenceBySite(outStartPosition)->getName());
      }
      outBotSegIt->setCoordinates(outStartPosition, inBotSegIt->getLength());
      for(hal_size_t inChild = 0; inChild < inNc; inChild++) {
        hal_size_t outChild = inChildToOutChild[inChild];
        if (outChild != outNc) {
          outBotSegIt->setChildIndex(outChild, inBotSegIt->getChildIndex(inChild));
          outBotSegIt->setChildReversed(outChild, inBotSegIt->getChildReversed(inChild));
        }
      }
      outBotSegIt->setTopParseIndex(inBotSegIt->getTopParseIndex());
    }
  }
}


void Genome::copySequence(Genome *dest) const
{
  DNAIteratorPtr inDna = getDNAIterator();
  DNAIteratorPtr outDna = dest->getDNAIterator();
  hal_size_t n = getSequenceLength();
  assert(n == dest->getSequenceLength());
  for (; (hal_size_t)inDna->getArrayIndex() < n; inDna->toRight(), 
         outDna->toRight())
  {
    outDna->setBase(inDna->getBase());
  }
  outDna->flush();
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
  BottomSegmentIteratorPtr botSegIt = getBottomSegmentIterator();
  TopSegmentIteratorPtr topSegIt = getTopSegmentIterator();
  int top = 0, bot = 0;
  while ((not botSegIt->atEnd()) && (not topSegIt->atEnd()))
  {
    bool bright = false;
    bool tright = false;
    BottomSegment* bseg = botSegIt->getBottomSegment();
    TopSegment* tseg = topSegIt->getTopSegment();
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
      botSegIt->toRight();
    }
    if (tright == true)
    {
      top += 1;
      topSegIt->toRight();
    }
  }
}
} // namespace hal
