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

namespace hal {

class LodNode;
class LodEdge;

class LodGraph
{
public:
   
   typedef std::multiset<LodNode*, LodNodePLess> NodeSet;
   typedef NodeSet::iterator NodeIterator;
   typedef std::map<const Sequence*, NodeSet*> NodeRefMap;
   typedef NodeRefMap::iterator NodeRefMapIterator;
   typedef std::map<const Genome*, NodeSet*> NodeMap;
   typedef NodeMap::iterator NodeMapIterator;

protected:

   
   // input alignment structure
   AlignmentConstPtr _alignment;
   const Genome* _parent;
   std::vector<const Genome*> _children;

   NodeRefMap _nodeRefMap;
   NodeMap _nodeMap;
};

}

#endif
