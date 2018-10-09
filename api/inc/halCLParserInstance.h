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

/** Get an instance of an HAL CLParser with the compile in storage
 * engine options added in.
 * @param createOptions set to true if tool will need to 
 * create a new hal file (as opposed to just reading an existing one)
 */
CLParserPtr halCLParserInstance(bool createOptions = false);
}

#endif
