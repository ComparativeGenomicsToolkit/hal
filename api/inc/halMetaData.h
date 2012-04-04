/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMETADATA_H
#define _HALMETADATA_H

#include <string>
#include <map>
#include "halDefs.h"

namespace hal {

/** 
 * Abstract base class for alignment (or genome) metadata
 * MetaData is a set of key/value pairs where each key and each 
 * value is represented by a string.
 */
class MetaData
{
public:
   virtual ~MetaData() = 0;
   
   virtual const std::map<std::string, std::string>& getMap() const = 0;
   virtual const std::string& get(const std::string& key) const = 0;
   virtual void set(const std::string& key, const std::string& value) = 0;
   virtual bool has(const std::string& key) const = 0;
};


}
#endif
