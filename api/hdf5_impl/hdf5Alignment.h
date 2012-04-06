/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5ALIGNMENT_H
#define _HDF5ALIGNMENT_H

#include <H5Cpp.h>
#include "halAlignment.h"
#include "halAlignmentInstance.h"
#include "halGenome.h"
#include "halMetaData.h"

namespace hal {

/** 
 * HDF5 implementation of hal::Alignment
 */
class HDF5Alignment : public Alignment
{
public:

   virtual ~HDF5Alignment();

   void createNew(const std::string& alignmentPath);
   void open(const std::string& alignmentPath, 
             bool readOnly);
   
   GenomePtr addGenome(const std::string& path, 
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
   // Nobody creates this class except through the interface. 
   friend AlignmentPtr hdf5AlignmentInstance();
   friend AlignmentConstPtr hdf5AlignmentInstanceReadOnly();
   friend AlignmentPtr 
   hdf5AlignmentInstance(const H5::FileCreatPropList& fileCreateProps,
                         const H5::FileAccPropList& fileAccessProps,
                         const H5::DSetCreatPropList& datasetCreateProps);
   friend AlignmentConstPtr 
   hdf5AlignmentInstanceReadOnly(const H5::FileCreatPropList& fileCreateProps,
                                 const H5::FileAccPropList& fileAccessProps,
                                 const H5::DSetCreatPropList& datasetCreateProps);

   HDF5Alignment();
   HDF5Alignment(const H5::FileCreatPropList& fileCreateProps,
                 const H5::FileAccPropList& fileAccessProps,
                 const H5::DSetCreatPropList& datasetCreateProps);

protected:

   H5::H5File* _file;
   H5::FileCreatPropList _cprops;
   H5::FileAccPropList _aprops;
   H5::DSetCreatPropList _dcprops;
   int _flags;
   MetaDataPtr _metaData;
   static const H5std_string MetaGroupName;
   
};

}
#endif

