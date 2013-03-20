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

namespace hal {

class LodEdge;
class LodGraph;
class LodNode;

/** Compare Lod Node pointers based on their startPositions */
struct LodNodePLess
{
   bool operator(const LodNode* n1, const LodNode* n2) const;
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
   typedef std::vector<LodEdge*> EdgeList;
   typedef EdgeList::iterator EdgeIterator;

public:

   LodNode();
   LodNode(const Sequence* sequence, hal_index_t start,
           hal_index_t last);
   ~LodNode();

   const Sequence* getSequence() const;
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;   
   
   /** Add a new edge between two nodes. Specifically from the
    * right side of src to the left side of target.   
    *
    * Implied constraint: src is left of target (on genome forward strand)
    */
   void addEdge(const Sequence* sequence, bool srcReversed, LodNode* tgt, 
                bool tgtReversed);

   /** Extend the node to by length bases.  If length is 0
    * then we extend as far as possible.  Reversed flag describes
    * which direction we extend (true = left, false = right).  */
   void extend(bool reversed, hal_size_t length = 0);
   
protected:

   /** Get the smallest edge arising from the given side side */
   EdgeIterator getMinEdge(bool reversed);

   const Sequence* _sequence;
   hal_index_t _startPosition;
   hal_index_t _endPostition;
   EdgeList _edges;

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
