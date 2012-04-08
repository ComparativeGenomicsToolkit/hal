/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTTEST_H
#define _HALALIGNMENTTEST_H

#include <vector>
#include "halAlignment.h"

struct TempCreateAlignment {
   TempCreateAlignment(hal::AlignmentPtr alignment);
   std::string _path;
   hal::AlignmentPtr _alignment;
   ~TempCreateAlignment();
};

struct TempReadAlignment {
   TempReadAlignment(hal::AlignmentPtr alignment,
                     const std::string& path);
   std::string _path;
   hal::AlignmentPtr _alignment;
   ~TempReadAlignment();
};

class AlignmentTest
{
   void check();   
   virtual void createCallBack(hal::AlignmentPtr alignment) = 0;
   virtual void checkCallBack(hal::AlignmentConstPtr alignment) = 0;
};

#endif
