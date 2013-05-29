/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBLOCKLIFTOVER_H
#define _HALBLOCKLIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "halLiftover.h"

namespace hal {

class BlockLiftover : public Liftover
{
public:
   
   BlockLiftover();
   virtual ~BlockLiftover();
                   
protected:

   void liftInterval();
   void visitBegin();

   void cleanTargetParalogies();
   
protected: 
   
   std::set<MappedSegmentConstPtr> _mappedSegments;
   SegmentIteratorConstPtr _refSeg;
   hal_index_t _lastIndex;
   std::set<const Genome*> _spanningTree;
};

}
#endif
