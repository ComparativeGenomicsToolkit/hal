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

   ~HDF5Alignment();

   void createNew(const std::string& alignmentPath);
   void open(const std::string& alignmentPath, 
             bool readOnly);
   void close();
   
   GenomePtr addGenome(const std::string& name,
                       const std::string& path,
                       const std::string& parentName,
                       const std::vector<std::string>& childNames);

   void removeGenome(const std::string& name);

   GenomeConstPtr openConstGenome(const std::string& name) const;

   GenomePtr openGenome(const std::string& name);

   std::string getRootName() const;

   std::string getParentName(const std::string& name) const;
   
   std::vector<std::string> 
   getChildNames(const std::string& name) const;

   std::vector<std::string>
   getLeafNamesBelow(const std::string& name) const;

   MetaDataPtr getMetaData();

   MetaDataConstPtr getMetaData() const;

   const std::string& getPath(const std::string& name) const;

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

