/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOMMON_H
#define _HALCOMMON_H

#include <map>
#include <string>
#include <vector>

namespace hal {

std::vector<std::string> chopString(const std::string& inString,
                                    const std::string& separator);

}
#endif

