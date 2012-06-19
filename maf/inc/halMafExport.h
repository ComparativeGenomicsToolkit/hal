/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFEXPORT_H
#define _HALMAFEXPORT_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"

namespace hal {

class MafExport
{
public:

   MafExport();
   void convert(hal::AlignmentConstPtr alignment,
                std::ostream& mafStream); 
   virtual ~MafExport();

protected:

   hal::AlignmentConstPtr _alignment;
   std::ostream* _mafStream;
};

}

#endif
