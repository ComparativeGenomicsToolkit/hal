/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5CLPARSER_H
#define _HDF5CLPARSER_H

#include <H5Cpp.h>
#include "halCLParser.h"

namespace hal {

/** 
 * Operations to specify HDF5 options in hal::CLParser and extract the
 * results.  Only provides function, not an object instance.
 */
class HDF5CLParser
{
public:
    static void defineOptions(CLParserPtr parser,
                              bool createOptions);
    static void applyToDCProps(CLParserPtr parser,
                               H5::DSetCreatPropList& dcprops);
    static void applyToAProps(CLParserPtr parser,
                              H5::FileAccPropList& aprops);
    static bool getInMemory(CLParserPtr parser);

    protected:
    friend class HDF5Alignment;

    private:
    // can't create
    HDF5CLParser() {
    }
};

}
#endif

