/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5ALIGNMENT_H
#define _HDF5ALIGNMENT_H

#include <H5Cpp.h>
#include "halAlignment.h"

namespace hal {

/** 
 * HDF5 implementation of hal::Alignment
 */
class HDF5Alignment : public Alignment
{
public:

   HDF5Alignment();
   virtual ~HDF5Alignment();

   void createNew(const std::string& alignmentPath);
   void open(const std::string& alignmentPath, 
                     bool readOnly);
   
   void addGenome(const std::string& path, 
                  const std::string* parentPath,
                  const std::vector<std::string>& childPaths);

   void removeGenome(const std::string& path);

   GenomeConstPtr openConstGenome(const std::string& path) const;

   GenomePtr openGenome(const std::string& path);

   std::string getParent(const std::string& path);
   
   std::vector<std::string> getChildren(const std::string& path);

   MetaDataPtr getMetaData();

   MetaDataConstPtr getMetaData() const;

protected:

   H5::H5File* _file;
   H5::FileCreatPropList _cprops;
   H5::FileAccPropList _aprops;
   int _flags;
   MetaDataPtr _metaData;
   static const H5std_string MetaGroupName;
   
};

}
#endif

