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

protected:

   typedef std::list<LodNode*> NodeList;
   typedef NodeList::iterator NodeIterator;
   typedef std::map<const Genome*, NodeList*> GenomeNodes;
   typedef GenomeNodes::iterator GenomeNodesIterator;

   void erase();
   void scanGenome(const Genome* genome, NodeList* nodeList);
protected:
   
   // input alignment structure
   AlignmentConstPtr _alignment;
   const Genome* _parent;
   std::vector<const Genome*> _children;

   // step size for interpolation
   hal_size_t _step;

   // ordered list of nodes for each sequence
   GenomeNodes _genomeNodes;

   // adjacency table used to generate the graph
   LodAdjTable _adjTable;
};

}

#endif
