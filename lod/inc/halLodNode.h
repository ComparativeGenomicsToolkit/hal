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
#include <map>
#include "hal.h"

namespace hal {

class LodEdge;
class LodGraph;
class LodNode;

/** Compare Lod Node pointers based on their startPositions */
struct LodNodePLess
{
   bool operator(const LodNode* n1, const LodNode* n2) const;
};

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
public:

   friend struct LodNodePLess;
   friend class LodGraph;
   typedef std::set<LodEdge*> EdgeSet;
   typedef EdgeSet::iterator EdgeIterator;

public:

   LodNode();
   LodNode(const Sequence* sequence, hal_index_t start,
           hal_index_t last);
   ~LodNode();

   const Sequence* getSequence() const;
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;   
   
   /** add a new edge between two nodes. specifically from the
    * right side of src to the left side of target.   */
   void addEdge(bool srcReversed, LodNode* tgt, bool tgtReversed);
protected:

   const Sequence* _sequence;
   hal_index_t _startPosition;
   hal_index_t _endPostition;
   EdgeSet _edgeSet;
   EdgeSet _edgeSetZL;

private:
   LodNode(const LodNode&);
   const LodNode& operator=(const LodNode&) const;
};

inline bool LodNodePLess::operator(const LodNode* n1, const LodNode* n2) const
{
  return n1->_startPosition < n2->_startPosition;
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


}

#endif
