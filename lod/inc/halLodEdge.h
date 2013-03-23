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

struct LodEdgePLess 
{
   bool operator()(const LodEdge* e1, const LodEdge* e2) const;
};

std::ostream& operator<<(std::ostream& os, const LodEdge& edge);

/** Edge class for the Level of Detail Graph.  Edges represent
 * gaps between the nodes along a particular sequence.  The goal
 * of the LOD graph is to reduce all edges to size 0, implying that
 * the nodes cover the original graph in its entirety.
 *
 * Edges can have different coordinates than their nodes.  These
 * coordinates are stored separately in the edge (pos1/2).  They
 * reflect the left and right endpoint of the edge.  The edge 
 * covers bases between these endpoints but not including them. 
 *
 * We enforce that _pos1 < _pos2
 */
class LodEdge
{
   friend std::ostream& operator<<(std::ostream& os, const LodEdge& edge);
   friend struct LodEdgePLess;
   friend class LodNode;
public:
   
   LodEdge();
   LodEdge(const Sequence* sequence, 
           LodNode* node1, hal_index_t pos1, bool reversed1,
           LodNode* node2, hal_index_t pos2, bool reversed2);
   ~LodEdge();

   const Sequence* getSequence() const;
   /* The start point of the edge relative to the node.  The edge
      contains bases (startPosition, endPosition) EXCLUSIVE */
   hal_index_t getStartPosition(const LodNode* node) const;

   /* The end point of the edge relative to the node.  The edge
      contains bases (startPosition, endPosition) EXCLUSIVE */
   hal_index_t getEndPosition(const LodNode* node) const;

   hal_size_t getLength() const;
   const LodNode* getNode1() const;
   const LodNode* getNode2() const;

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

   void shrink(hal_size_t delta, bool fromLeft);

   const Sequence* _sequence;

   // genome coordinates
   LodNode* _node1;
   hal_index_t _pos1;
   LodNode* _node2;
   hal_index_t _pos2;
   bool _reversed1;
   bool _reversed2;
};

inline const Sequence* LodEdge::getSequence() const
{
  return _sequence;
}

inline hal_index_t LodEdge::getStartPosition(const LodNode* node) const
{
  assert(node == _node1 || node == _node2);
  return node == _node1 ? _pos1 : _pos2;
}
inline hal_index_t LodEdge::getEndPosition(const LodNode* node) const
{
  assert(node == _node1 || node == _node2);
  return node == _node1 ? _pos2 : _pos1;
}

inline hal_size_t LodEdge::getLength() const
{
  return _pos2 - _pos1 - 1;
}

inline const LodNode* LodEdge::getNode1() const
{
  return _node1;
}

inline const LodNode* LodEdge::getNode2() const
{
  return _node2;
}

}

#endif
