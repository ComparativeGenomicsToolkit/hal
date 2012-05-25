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
  cout << _tree << endl << endl;

  cout << "GenomeName, NumChildren, Length, NumSequences, "
       << "NumTopSegments, NumBottomSegments" << endl;

  vector<GenomeStats>::const_iterator i;
  for (i = _genomeStatsVec.begin(); i != _genomeStatsVec.end(); ++i)
  {
    cout << i->_name << ", " << i->_numChildren << ", " << i->_length
         << ", " << i->_numSequences << ", " << i->_numTopSegments
         << ", " << i->_numBottomSegments << endl;
  }
  cout << endl;
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
    readGenomeRecursive(root);
  }
}

void HalStats::readGenomeRecursive(const Genome* genome)
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

  for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
  {
    const Genome* child = genome->getChild(i);
    readGenomeRecursive(child);
  }
}
