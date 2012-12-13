/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMASKEXTRACTOR_H
#define _HALMASKEXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"

namespace hal {

class MaskExtractor 
{
public:

   MaskExtractor();
   virtual ~MaskExtractor();

   void extract(AlignmentConstPtr alignment, std::ostream* bedStream, 
                hal_size_t extend, double extendPct);

protected:

   AlignmentConstPtr _alignment;
   std::ostream* _bedStream;
   hal_size_t _extend; 
   double _extendPct;
   PositionCache _posCache;
   
};

}


#endif
