/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCLPARSERINSTANCE_H
#define _HALCLPARSERINSTANCE_H

#include "halDefs.h"
#include "halCLParser.h"

/** HAL API namespace */
namespace hal {

/** Get an instance of an HDF5-implemented CLParser 
 * @param createOptions set to true if tool will need to 
 * create a new hal file (as opposed to just reading an existing one)
 */
CLParserPtr hdf5CLParserInstance(bool createOptions = false);

// if any other implementations are added, could add functionality
// to detect which types of instance to add from command line options
// (ie --hdf5 or something) but not gonna waste any time on that now.
}

#endif
