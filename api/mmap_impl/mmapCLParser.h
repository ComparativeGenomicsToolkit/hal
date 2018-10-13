
/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _MMAPCLPARSER_H
#define _MMAPCLPARSER_H

#include "halCLParser.h"

namespace hal {

/** 
 * Operations to specify mmap options in hal::CLParser and extract the
 * results.  Only provides function, not an object instance.
 */
class MMapCLParser
{
public:
    static void defineOptions(CLParserPtr parser,
                              unsigned mode);
    static size_t getInitSize(CLParserPtr parser);
    static size_t getGrowSize(CLParserPtr parser);

    
    protected:
   //friend class MMapAlignment;

    private:
    // can't create
    MMapCLParser() {
    }
};

}
#endif

