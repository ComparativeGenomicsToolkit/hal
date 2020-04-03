/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBEDSCANNER_H
#define _HALBEDSCANNER_H

#include "hal.h"
#include "halBedLine.h"
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace hal {

    /** Parse a BED file line by line
     * written independently from the bed export, and it's too much of a
     * bother to reuse any of that code. */
    class BedScanner {
      public:
        BedScanner();
        virtual ~BedScanner();
        virtual void scan(const std::string &bedPath, int bedType=0);
        virtual void scan(std::istream *bedStream, int bedType=0);

        static size_t getNumColumns(const std::string &bedLine);

      protected:
        virtual void visitBegin();
        virtual void visitLine();
        virtual void visitEOF();

        static void skipWhiteSpaces(std::istream *bedStream);

      protected:
        std::istream *_bedStream;
        BedLine _bedLine;
        hal_size_t _lineNumber;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
