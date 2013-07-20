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
    Sequence::UpdateInfo info(sequence->getName(),
                              sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
  dest->updateBottomDimensions(dimensions);
}

void Genome::copyTopSegments(Genome *dest) const
{
  TopSegmentIteratorConstPtr inTop = getTopSegmentIterator();
  TopSegmentIteratorPtr outTop = dest->getTopSegmentIterator();
  hal_size_t n = dest->getNumTopSegments();
  assert(n == 0 || n == getNumTopSegments());
  for (; (hal_size_t)inTop->getArrayIndex() < n; inTop->toRight(),
         outTop->toRight())
  {
    outTop->setCoordinates(inTop->getStartPosition(), inTop->getLength());
    outTop->setParentIndex(inTop->getParentIndex());
    outTop->setParentReversed(inTop->getParentReversed());
    outTop->setBottomParseIndex(inTop->getBottomParseIndex());
    outTop->setNextParalogyIndex(inTop->getNextParalogyIndex());
  }
}

void Genome::copyBottomSegments(Genome *dest) const
{
  BottomSegmentIteratorPtr inBot = getBottomSegmentIterator();
  BottomSegmentIteratorPtr outBot = dest->getBottomSegmentIterator();
  hal_size_t n = getNumBottomSegments();
  hal_size_t inNc = getNumChildren();
  hal_size_t outNc = dest->getNumChildren();
  assert(n == dest->getNumBottomSegments());
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
  for (; (hal_size_t)inBot->getArrayIndex() < n; inBot->toRight(),
         outBot->toRight())
  {
    outBot->setCoordinates(inBot->getStartPosition(), inBot->getLength());
    for(hal_size_t inChild = 0; inChild < inNc; inChild++) {
      hal_size_t outChild = inChildToOutChild[inChild];
      if(outChild != outNc) {
        outBot->setChildIndex(outChild, inBot->getChildIndex(inChild));
        outBot->setChildReversed(outChild, inBot->getChildReversed(inChild));
      }
    }
    outBot->setTopParseIndex(inBot->getTopParseIndex());
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
