/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODEXTRACT_H
#define _HALLODEXTRACT_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"
#include "halLodGraph.h"

namespace hal {


class LodExtract
{
public:
   
   LodExtract();
   ~LodExtract();

   void createInterpolatedAlignment(AlignmentConstPtr inAlignment,
                                    AlignmentPtr outAlignment,
                                    hal_size_t step,
                                    const std::string& tree);
   
   
protected:

   void createTree(const std::string& tree);
   void convertInternalNode(const std::string& genomeName);
   

   AlignmentConstPtr _inAlignment;
   AlignmentPtr _outAlignment;

   LodGraph _graph;
};

}

#endif
