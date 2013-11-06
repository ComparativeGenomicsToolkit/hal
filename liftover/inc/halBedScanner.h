/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBEDSCANNER_H
#define _HALBEDSCANNER_H

#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <locale>
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
   virtual void scan(const std::string& bedPath, int bedVersion = -1,
                     const std::locale* inLocale = NULL);
   virtual void scan(std::istream* bedStream, int bedVersion = -1,
                     const std::locale* inLocale = NULL);
   
   static int getBedVersion(std::istream* bedStream, 
                            const std::locale* inLocale = NULL);
   static size_t getNumColumns(const std::string& bedLine,
                               const std::locale* inLocale = NULL);

protected:
   
   virtual void visitBegin();
   virtual void visitLine();
   virtual void visitEOF();
   
   static void skipWhiteSpaces(std::istream* bedStream,
                               const std::locale* inLocale = NULL);

protected:

   std::istream* _bedStream;
   BedLine _bedLine;
   hal_size_t _lineNumber;
   int _bedVersion;
};

}

#endif
