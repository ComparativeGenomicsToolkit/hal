/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENT_H
#define _HALALIGNMENT_H

#include <string>
#include "halDefs.h"
#include "halMetaData.h"
#include "halGenome.h"

namespace hal {

/** 
 * Abstract base class for a hierarhcical alignment
 */
class Alignment
{
public:

   Alignment();
   virtual ~Alignment();

   virtual void createNew(const std::string& alignmentPath) = 0;
   virtual void open(const std::string& alignmentPath, 
                     bool readOnly) = 0;
   
   virtual void addGenome(const std::string& path, 
                          const std::string* parentPath,
                          const std::vector<std::string>& childPaths) = 0;

   virtual void removeGenome(const std::string& path) = 0;

   virtual GenomeConstPtr openConstGenome(const std::string& path) const = 0;

   virtual GenomePtr openGenome(const std::string& path) = 0;

   virtual std::string getParent(const std::string& path) = 0;
   
   virtual std::vector<std::string> getChildren(const std::string& path) = 0;

   virtual MetaDataPtr getMetaData() = 0;

   virtual MetaDataConstPtr getMetaData() const = 0;
};


}
#endif
