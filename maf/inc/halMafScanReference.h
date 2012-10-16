/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFSCANREFERENCE_H
#define _HALMAFSCANREFERENCE_H

#include <iostream>
#include <string>
#include <deque>
#include <cstdlib>
#include <map>
#include <string>
#include "halMafScanner.h"

namespace hal {

/** Parse a MAF file line by line, getting some dimension stats
 * and maybe checking for some errros. */
class MafScanReference : private MafScanner
{
public:

   MafScanReference();
   ~MafScanReference();
   
   std::string getRefName(const std::string& mafPath);
   
private:

   void aLine();
   void sLine();
   void end();

   std::string _name;
};

}

#endif
