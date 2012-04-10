/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5ALIGNMENT_H
#define _HDF5ALIGNMENT_H

#include <map>
#include <H5Cpp.h>
#include "halAlignment.h"
#include "halAlignmentInstance.h"
#include "halGenome.h"
#include "halMetaData.h"

typedef struct _stTree stTree;

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
                       const std::pair<std::string, double>& parentName,
                       const std::vector<std::pair<std::string, double> >&
                       childNames);

   void removeGenome(const std::string& name);

   GenomeConstPtr openConstGenome(const std::string& name) const;

   GenomePtr openGenome(const std::string& name);

   std::string getRootName() const;

   std::string getParentName(const std::string& name) const;

   double getBranchLength(const std::string& parentName,
                          const std::string& childName);   

   std::vector<std::string> 
   getChildNames(const std::string& name) const;

   std::vector<std::string>
   getLeafNamesBelow(const std::string& name) const;

   hal_size_t getNumGenomes() const;

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
   
   void loadTree();
   void writeTree();
   void addGenomeToTree(const std::string& name,
                        const std::pair<std::string, double>& parentName,
                        const std::vector<std::pair<std::string, double> >&
                        childNames);

protected:

   H5::H5File* _file;
   H5::FileCreatPropList _cprops;
   H5::FileAccPropList _aprops;
   H5::DSetCreatPropList _dcprops;
   int _flags;
   MetaDataPtr _metaData;
   static const H5std_string MetaGroupName;
   static const H5std_string TreeGroupName;
   static const H5std_string GenomesGroupName;
   stTree* _tree;
   bool _dirty;
   mutable std::map<std::string, GenomePtr> _openGenomes;
};

}
#endif

