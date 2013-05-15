/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBEDSCANNER_H
#define _HALBEDSCANNER_H

#include <fstream>
#include <string>
#include <deque>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include "hal.h"
#include "halBedLine.h"

namespace hal {

/** Parse a BED file line by line 
 * written independently from the bed export, and it's too much of a 
 * bother to reuse any of that code. */
class BedScanner
{
public:
   BedScanner();
   virtual ~BedScanner();
   virtual void scan(const std::string& bedPath, int bedVersion = -1);
   
   static int getBedVersion(const std::string& bedPath);

protected:

   virtual void visitLine();
   virtual void visitEOF();

protected:

   std::ifstream _bedFile;
   BedLine _bedLine;
};

}

#endif
