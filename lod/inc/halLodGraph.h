/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODGRAPH_H
#define _HALLODGRAPH_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"
#include "halLodSegment.h"
#include "halLodBlock.h"

namespace hal {

class LodGraph
{
public:

   typedef std::set<LodSegment*, LodSegmentPLess> SegmentSet;
   typedef SegmentSet::iterator SegmentIterator;
   
   LodGraph();
   ~LodGraph();

   void erase();

   const LodBlock* getBlock(hal_size_t index) const;
   hal_size_t getNumBlocks() const;
   const SegmentSet* getSegmentSet(const Sequence* sequence) const;
   const LodBlock* getTelomeres() const;

   /** Build the LOD graph for a given subtree of the alignment.  The
    * entire graph is stored in memory in a special structure (ie not within
    * HAL).  The step parameter dictates how coarse-grained the interpolation
    * is:  every step bases are sampled.  */
   void build(AlignmentConstPtr alignment, const Genome* parent,
              const std::vector<const Genome*>& children, 
              const Genome* grandParent,
              hal_size_t step, bool allSequences, double probeFrac,
              double minSeqFrac);

   /** Help debuggin and tuning */
   void printDimensions(std::ostream& os) const;

   /** Make sure that every base in every genome appears exactly once
    * in the sequence graph.  For debugging purpose only -- will take
    * up a significant amount of memory and time */
   bool checkCoverage() const;

protected:

   typedef std::vector<LodBlock*> BlockList;
   typedef BlockList::iterator BlockIterator;
   typedef BlockList::const_iterator BlockConstIterator;

   typedef std::map<const Sequence*, SegmentSet*> SequenceMap;
   typedef SequenceMap::iterator SequenceMapIterator;

   /** Read a HAL genome into sequence graph */
   void scanGenome(const Genome* genome);

   /** Check maxium distance of this column to any other sampled position.
    * Also count the number of genomes it aligns to.  This information
    * will be used to prioritize probed columns*/
   void evaluateColumn(ColumnIteratorConstPtr colIt, hal_size_t& outDeltaMax,
                       hal_size_t& outNumGenomes, hal_size_t& outMinSeqLen);

   /* Test if this is the best column based on stats collected above */
   bool bestColumn(hal_size_t probeStep, hal_size_t delta, 
                   hal_size_t numGenomes, hal_size_t minSeqLen,
                   hal_size_t maxDelta, hal_size_t maxNumGenomes,
                   hal_size_t maxMinSeqLen);

   /** Add segments for the telomeres of the sequence, ie at 
    * position -1 and and endPosition + 1 */
   void addTelomeres(const Sequence* sequence);

   /** Add a single column iterator as  a block */
   void createColumn(ColumnIteratorConstPtr colIt);

   /** Add an entire sequence as unaliged segment */
   void createUnaligedSegment(const Sequence* sequence);

   /** compute the adjacencies using the SegmentSets */
   void computeAdjacencies();

   /** First optimization pass: Maximally extend all blocks */
   void optimizeByExtension();

   /** Second optimization pass: Merge all compatible adjacent blocks */
   void optimizeByMerging();

   /** Third optimization pass: Insert new blocks until all edges have 
    * zero length */
   void optimizeByInsertion();

protected:
   
   // input alignment structure
   AlignmentConstPtr _alignment;
   const Genome* _parent;
   std::set<const Genome*> _genomes;
   const Genome* _grandParent;

   // step size for interpolation
   hal_size_t _step;

   // fraction of edge to greedily extend
   double _extendFraction;

   // the alignment blocks
   BlockList _blocks;

   // the telomeres all get put in one block.  the block structure
   // is used only to make sure they get freed.
   LodBlock _telomeres;
   
   // nodes sorted by sequence
   SequenceMap _seqMap;

   // sample all sequences no matter how small they are
   bool _allSequences;

   // number of probe positions to try in range as fraction of the
   // step size
   // [pos - step / 2, pos + step / 2] before accepting column
   double _probeFrac;
   
   // min size of sequence to not be ignored (computed from the
   // minSeqFrac paramater)
   hal_size_t _minSeqLen;
};

inline const LodBlock* LodGraph::getBlock(hal_size_t index) const
{
  return _blocks[index];
}

inline hal_size_t LodGraph::getNumBlocks() const
{
  return _blocks.size();
}

inline const LodGraph::SegmentSet* LodGraph::getSegmentSet(
  const Sequence* sequence) const
{
  assert(_seqMap.find(sequence) != _seqMap.end());
  return _seqMap.find(sequence)->second;
}

inline const LodBlock* LodGraph::getTelomeres() const
{
  return &_telomeres;
}

}

#endif
