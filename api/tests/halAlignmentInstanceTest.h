/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTINSTANCETEST_H
#define _HALALIGNMENTINSTANCETEST_H

#include "halAlignmentInstance.h"
#include <vector>

using namespace hal;
Alignment *getTestAlignmentInstances(const std::string &storageFormat, const std::string &alignmentPath, unsigned mode);

#endif
// Local Variables:
// mode: c++
// End:
