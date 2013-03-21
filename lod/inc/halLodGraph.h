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
#include "halLodAdjTable.h"

namespace hal {

class LodNode;

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
              const std::vector<const Genome*> children, 
              hal_size_t step);

   /** Help debuggin and tuning */
   void printDimensions(std::ostream& os) const;

protected:

   typedef std::list<LodNode*> NodeList;
   typedef NodeList::iterator NodeIterator;
   typedef std::map<const Genome*, NodeList*> GenomeNodes;
   typedef GenomeNodes::iterator GenomeNodesIterator;

   void erase();

   /** Read a HAL genome into a nodeList, updating the adjTable as well */
   void scanGenome(const Genome* genome, NodeList* nodeList);

   /** Add a single column iterator */
   void createColumn(ColumnIteratorConstPtr colIt, NodeList* nodeList);

   /** First optimization pass: Maximally extend all nodes */
   void optimizeByExtension();

   /** Second optimization pass: Insert new nodes until all edges have 
    * zero length */
   void optimizeByInsertion();

protected:

   
   // input alignment structure
   AlignmentConstPtr _alignment;
   const Genome* _parent;
   std::vector<const Genome*> _children;

   // step size for interpolation
   hal_size_t _step;

   // fraction of edge to greedily extend
   double _extendFraction;

   // ordered list of nodes for each sequence
   GenomeNodes _genomeNodes;

   // adjacency table used to generate the graph
   LodAdjTable _adjTable;
};

}

#endif
