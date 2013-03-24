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
   
   LodGraph();
   ~LodGraph();

   /** Build the LOD graph for a given subtree of the alignment.  The
    * entire graph is stored in memory in a special structure (ie not within
    * HAL).  The step parameter dictates how coarse-grained the interpolation
    * is:  every step bases are sampled.  */
   void build(AlignmentConstPtr alignment, const Genome* parent,
              const std::vector<const Genome*>& children, 
              hal_size_t step);

   /** Help debuggin and tuning */
   void printDimensions(std::ostream& os) const;

protected:

   typedef std::vector<LodBlock*> BlockList;
   typedef BlockList::iterator BlockIterator;
   typedef BlockList::const_iterator BlockConstIterator;

   typedef std::set<LodSegment*, LodSegmentPLess> SegmentSet;
   typedef SegmentSet::iterator SegmentIterator;

   typedef std::map<const Sequence*, SegmentSet*> SequenceMap;
   typedef SequenceMap::iterator SequenceMapIterator;

   void erase();

   /** Read a HAL genome into sequence graph */
   void scanGenome(const Genome* genome);

   /** Test if we can add a column.  Does it collide?  does it fail 
    * heuristics? */
   bool canAddColumn(ColumnIteratorConstPtr colIt);

   /** Add segments for the telomeres of the sequence, ie at 
    * position -1 and and endPosition + 1 */
   void addTelomeres(const Sequence* sequence);

   /** Add a single column iterator as  a block */
   void createColumn(ColumnIteratorConstPtr colIt);

   /** compute the adjacencies using the SegmentSets */
   void computeAdjacencies();

   /** First optimization pass: Maximally extend all blocks */
   void optimizeByExtension();

   /** Second optimization pass: Insert new blocks until all edges have 
    * zero length */
   void optimizeByInsertion();

protected:

   
   // input alignment structure
   AlignmentConstPtr _alignment;
   const Genome* _parent;
   std::set<const Genome*> _genomes;

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
};

}

#endif
