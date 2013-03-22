/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODNODE_H
#define _HALLODNODE_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <map>
#include "hal.h"
#include "halLodEdge.h"

namespace hal {

class LodGraph;
class LodNode;

/** Compare Lod Node pointers based on their startPositions */
struct LodNodePLess
{
   bool operator()(const LodNode* n1, const LodNode* n2) const;
};

std::ostream& operator<<(std::ostream& os, const LodNode& node);

/* Node strcture for the Level of Detail graph.  These Nodes are what 
 * will eventually be turned into (top or bottom) segments in the output
 * HAL graph.  Nodes have two edge sets (bidirected graph) one for the 
 * head and one for the tail.  Each node is associated with a segment 
 * on a particular genome, but its adjacency list is governed by its 
 * homologous positions on all other genomes in consideration.  
 *
 * Note that _endPosition is the last postion (ie not last + 1)
 *
 * All edges incident to the node get freed by the destructor.
 */
class LodNode
{
   friend std::ostream& operator<<(std::ostream& os, const LodNode& node);

public:

   friend struct LodNodePLess;
   friend class LodGraph;

public:

   LodNode();
   LodNode(const Sequence* sequence, hal_index_t start,
           hal_index_t last);
   ~LodNode();

   const Sequence* getSequence() const;
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;
   hal_size_t getDegree() const;
   
   /** Add a new edge between two nodes. Specifically from the
    * right side of src to the left side of target.   
    *
    * Implied constraint: src is left of target (on genome forward strand)
    */
   void addEdge(const Sequence* sequence,
                hal_index_t srcPos, bool srcReversed, 
                LodNode* tgt, hal_index_t tgtPos, bool tgtReversed);

   /** Extend the node by extendFraction of the maximum possible 
    * number of bases in both directions, starting with forward */
   void extend(double extendFraction = 1.);

   /** Remove edge length by creating a new node for each edge
    * and storing them in the buffer.  Should call after extend */
   void fillInEdges(std::vector<LodNode*>& nodeBuffer);   
   
   /** Get edge length stats.  Most useful for debugging.  If just
    * need min lengths then getMinLengths is better.*/
   void getEdgeLengthStats(hal_size_t& fMin, hal_size_t& fMax, 
                           hal_size_t& fTot,
                           hal_size_t& rMin, hal_size_t& rMax,
                           hal_size_t& rTot) const;
   
protected:

   typedef std::vector<LodEdge*> EdgeList;
   typedef EdgeList::iterator EdgeIterator;
   typedef EdgeList::const_iterator EdgeConstIterator;

   /* Find the minimum edge lengths in both directions */
   void getMinLengths(hal_size_t& rMin, hal_size_t& fMin) const;

   /** Insert a node into all edges in the batch.  These 
    * edges must be compatible in terms of their sequence locations. 
    * ie each edge must represent the exact same stretch of sequence */
   void insertFillNode(std::vector<LodEdge*>& edgeBatch,
                       std::vector<LodNode*>& nodeBuffer);

   const Sequence* _sequence;
   hal_index_t _startPosition;
   hal_index_t _endPosition;
   EdgeList _edges;

private:
   LodNode(const LodNode&);
   const LodNode& operator=(const LodNode&) const;
};

inline bool LodNodePLess::operator()(const LodNode* n1, const LodNode* n2) const
{
  return n1->_startPosition < n2->_startPosition;
}

inline const Sequence* LodNode::getSequence() const
{
  return _sequence;
}

inline hal_index_t LodNode::getStartPosition() const 
{
  return _startPosition;
}

inline hal_index_t LodNode::getEndPosition() const
{
  return _endPosition;
}

inline hal_size_t LodNode::getLength() const
{
  hal_index_t length = 1 + _endPosition - _startPosition;
  assert (length > 0);
  return (hal_size_t)length;
}

inline hal_size_t LodNode::getDegree() const
{
  return _edges.size();
}

}

#endif
