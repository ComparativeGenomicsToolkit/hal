/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODADJTABLE_H
#define _HALLODADJTABLE_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"

namespace hal {

class LodEdge;
class LodNode;

/** This is a temporary structure used to scan in all the adjacencies that
 * are used to construct the LodGraph.  The final adjacencies are then
 * written back into the nodes.  This structure could be preserved throughout
 * but am worried that accesses would be too slow. 
 *
 * This structure contains no invented information.  It exactly reflects
 * the original graph properties (ie no interpolating is done).  
 *
 * Nodes cannot overlap in this structure.  If some node has been added that
 * covers position X in sequence Y, any new node added that covers this position
 * is skipped.  
 *
 * All coordinates are in HAL internal format (forward genome)
 *
 * It's critical to the current implementation that all sequence endpoints
 * get added to the graph.
 *
 * For each seqeunce, store a sorted list of "Node References" which are 
 * nodes paired with the positions they are homolgous to on that sequence
 */
class LodAdjTable
{
public:
   LodAdjTable();
   ~LodAdjTable();
   
   /** Add a new node (which has already been created) to the 
    * adjacency table, using the column iterator to resolve all
    * its homologous positions.  */
   void addNode(LodNode* node, ColumnIteratorConstPtr colIt);

   /** For every sequence, for every pair of adjacent node refs,
    * add an edge connecting the two nodes */
   void writeAdjacenciesIntoNodes();

   /** Clear all internally created structures.  Note that all LodNode*
    * elements (added with addNode) are of course not freed.  Nor are
    * any edges that were added to them */
   void clear();

protected:
   
   /** An instance of a node positioned along a sequence */
   struct NodeRef 
   {
      NodeRef(hal_index_t pos = 0, bool reversed = false, 
              LodNode* _node = NULL);
      hal_index_t _pos;
      bool _reversed;
      LodNode* _node;
      bool operator<(const NodeRef& other) const;
   };
   
   typedef std::multiset<NodeRef> RefSet;
   typedef RefSet::iterator RefIterator;
   typedef std::map<const Sequence*, RefSet*> Table;
   typedef Table::iterator Iterator;

protected:

   /** the adjacency table */
   Table _table;
};

inline LodAdjTable::NodeRef::NodeRef(hal_index_t pos, bool reversed, 
                                     LodNode* node) : 
  _pos(pos), _reversed(reversed), _node(node)
{
}

inline bool LodAdjTable::NodeRef::operator<(const NodeRef& other) const
{
  return _pos < other._pos;
}

}

#endif
