#include "hal.h"

using namespace std;
using namespace hal;

// copied from submodules/hal/extract/halExtract.cpp
void copyAllDimensions(AlignmentConstPtr inAlignment, const Genome* inGenome,
                       Genome *outGenome)
{
  vector<Sequence::Info> dimensions;
  SequenceIteratorConstPtr seqIt = inGenome->getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = inGenome->getSequenceEndIterator();

  bool root = inAlignment->getParentName(inGenome->getName()).empty();
  bool leaf = inAlignment->getChildNames(inGenome->getName()).empty();     
  
  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    Sequence::Info info(sequence->getName(),
                        sequence->getSequenceLength(),
                        root ? 0 : sequence->getNumTopSegments(),
                        leaf ? 0 : sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
  outGenome->setDimensions(dimensions);
}

// mostly copied from submodules/hal/extract/halExtract.cpp
void copyBotDimensions(const Genome *inGenome, Genome *outGenome)
{
  vector<Sequence::UpdateInfo> dimensions;
  SequenceIteratorConstPtr seqIt = inGenome->getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = inGenome->getSequenceEndIterator();

  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    Sequence::UpdateInfo info(sequence->getName(),
                              sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
  outGenome->updateBottomDimensions(dimensions);
}

void copyTopDimensions(const Genome *inGenome, Genome *outGenome)
{
  vector<Sequence::UpdateInfo> dimensions;
  SequenceIteratorConstPtr seqIt = inGenome->getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = inGenome->getSequenceEndIterator();
  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    Sequence::UpdateInfo info(sequence->getName(),
                              sequence->getNumTopSegments());
    dimensions.push_back(info);
  }
  outGenome->updateTopDimensions(dimensions);
}

// mostly copied from submodules/hal/extract/halExtract.cpp
void copyGenomeWithoutBotSegments(const Genome *inGenome, Genome *outGenome)
{
  DNAIteratorConstPtr inDna = inGenome->getDNAIterator();
  DNAIteratorPtr outDna = outGenome->getDNAIterator();
  hal_size_t n = inGenome->getSequenceLength();
  assert(n == outGenome->getSequenceLength());
  for (; (hal_size_t)inDna->getArrayIndex() < n; inDna->toRight(), 
         outDna->toRight())
  {
    outDna->setChar(inDna->getChar());
  }

  TopSegmentIteratorConstPtr inTop = inGenome->getTopSegmentIterator();
  TopSegmentIteratorPtr outTop = outGenome->getTopSegmentIterator();
  n = outGenome->getNumTopSegments();
  assert(n == 0 || n == inGenome->getNumTopSegments());
  for (; (hal_size_t)inTop->getArrayIndex() < n; inTop->toRight(),
         outTop->toRight())
  {
    outTop->setCoordinates(inTop->getStartPosition(), inTop->getLength());
    outTop->setParentIndex(inTop->getParentIndex());
    outTop->setParentReversed(inTop->getParentReversed());
    outTop->setBottomParseIndex(inTop->getBottomParseIndex());
    outTop->setNextParalogyIndex(inTop->getNextParalogyIndex());
  }

  const map<string, string>& meta = inGenome->getMetaData()->getMap();
  map<string, string>::const_iterator i = meta.begin();
  for (; i != meta.end(); ++i)
  {
    outGenome->getMetaData()->set(i->first, i->second);
  }
}

void copyBotSegments(const Genome *inGenome, Genome *outGenome)
{
  BottomSegmentIteratorPtr inBot = inGenome->getBottomSegmentIterator();
  BottomSegmentIteratorPtr outBot = outGenome->getBottomSegmentIterator();
  hal_size_t n = inGenome->getNumBottomSegments();
  hal_size_t inNc = inGenome->getNumChildren();
  hal_size_t outNc = outGenome->getNumChildren();

  // The child indices aren't consistent across files--make sure each bottom
  // segment points to the correct children
  vector<string> inChildNames;
  vector<string> outChildNames;
  for (hal_size_t inChild = 0; inChild < inNc; ++inChild)
  {
    inChildNames.push_back(inGenome->getChild(inChild)->getName());
  }
  for (hal_size_t outChild = 0; outChild < outNc; ++outChild)
  {
    outChildNames.push_back(outGenome->getChild(outChild)->getName());
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

void fixParseInfo(Genome *genome)
{
  if (genome->getParent() == NULL || genome->getNumChildren() == 0)
  {
    return;
  }
  
  // copied from CactusHalConverter::updateRootParseInfo() in
  // cactus2hal/src/cactusHalConverter.cpp 
  BottomSegmentIteratorPtr bottomIterator = 
    genome->getBottomSegmentIterator();
  TopSegmentIteratorPtr topIterator = genome->getTopSegmentIterator();
  BottomSegmentIteratorConstPtr bend = genome->getBottomSegmentEndIterator();
  TopSegmentIteratorConstPtr tend = genome->getTopSegmentEndIterator();
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
