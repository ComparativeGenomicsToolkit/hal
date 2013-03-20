/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODEDGE_H
#define _HALLODEDGE_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"

namespace hal {

class LodNode;
class LodEdge;

std::ostream& operator<<(std::ostream& os, const LodEdge& edge);

/** Edge class for the Level of Detail Graph.  Edges represent
 * gaps between the nodes along a particular sequence.  The goal
 * of the LOD graph is to reduce all edges to size 0, implying that
 * the nodes cover the original graph in its entirety 
 *
 * Nodes in the edge are not stored arbitrarily:  Node1 is has the 
 * address with the smaller value.  
 */
class LodEdge
{
   friend std::ostream& operator<<(std::ostream& os, const LodEdge& edge);

public:
   
   LodEdge();
   LodEdge(const Sequnce* sequence, size_t length, LodNode* node1,
           bool reversed1, LodNode* node2, bool reversed2);
   ~LodEdge();
   
   const Sequence* getSeqeunce() const;
   hal_size_t getLength() const;

   /** Given one node of the edge, get the other node. 
    * @param node The input node (must be in the edge)
    * @param revThis If not null, then orientation of this is stored here
    * @param revOther If not null, then orientation of other stored here
    * @return the other node
    */
   LodNode* getOtherNode(const LodNode* node, bool* revThis = NULL,
                         bool* revOther = NULL);

   /** Replace input node pointer with NULL in the edge.  This is used
    * only as a hack to speed up the node destructor */
   void nullifyNode(const LodNode* node);

   

protected:

   const Sequencey* _sequence;
   hal_size_t _length;

   bool _reversed1;
   LodNode* _node1;
   bool _reversed2;
   LodNode* _node2;   
};



inline const Sequence* LodEdge::getSequence() const
{
  return _sequence;
}

inline hal_size_t LodEdge::getLength() const
{
  return _length;
}

inline bool LodEdge::getOtherNode(const LodNode* node) const
{
  assert(node == _node1 || node == _node2);
  return node == _node
}

}

#endif
