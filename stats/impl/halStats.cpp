/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halStats.h"

using namespace std;
using namespace hal;

HalStats::HalStats()
{

}

HalStats::HalStats(AlignmentConstPtr alignment)
{
  readAlignment(alignment);
}

HalStats::~HalStats()
{

}

void HalStats::printCsv(ostream& outStream) const
{
  outStream << _tree << endl << endl;

  outStream << "GenomeName, NumChildren, Length, NumSequences, "
       << "NumTopSegments, NumBottomSegments" << endl;

  vector<GenomeStats>::const_iterator i;
  for (i = _genomeStatsVec.begin(); i != _genomeStatsVec.end(); ++i)
  {
    outStream << i->_name << ", " << i->_numChildren << ", " << i->_length
         << ", " << i->_numSequences << ", " << i->_numTopSegments
         << ", " << i->_numBottomSegments << endl;
  }
  outStream << endl;
}

void HalStats::readAlignment(AlignmentConstPtr alignment)
{
  _tree.clear();
  _genomeStatsVec.clear();

  if (alignment->getNumGenomes() > 0)
  {
    _tree = alignment->getNewickTree();
    _genomeStatsVec.reserve(alignment->getNumGenomes());
    const Genome* root = alignment->openGenome(alignment->getRootName());
    readGenomeRecursive(alignment, root);
  }
}

void HalStats::readGenomeRecursive(AlignmentConstPtr alignment,
                                   const Genome* genome)
{
  assert(genome != NULL);

  GenomeStats genomeStats;
  genomeStats._name = genome->getName();
  genomeStats._numChildren = genome->getNumChildren();
  genomeStats._length = genome->getSequenceLength();
  genomeStats._numSequences = genome->getNumSequences();
  genomeStats._numTopSegments = genome->getNumTopSegments();
  genomeStats._numBottomSegments = genome->getNumBottomSegments();
  _genomeStatsVec.push_back(genomeStats);

  vector<string> children = alignment->getChildNames(genome->getName());
  for (hal_size_t i = 0; i < children.size(); ++i)
  {
    const Genome* child = alignment->openGenome(children[i]);
    readGenomeRecursive(alignment, child);
  }
}
