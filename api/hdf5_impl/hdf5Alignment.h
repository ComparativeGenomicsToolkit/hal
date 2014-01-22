/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5ALIGNMENT_H
#define _HDF5ALIGNMENT_H

#include <map>
#include <H5Cpp.h>
#include "hdf5Alignment.h"
#include "halAlignmentInstance.h"
#include "hdf5Genome.h"
#include "hdf5MetaData.h"

typedef struct _stTree stTree;

namespace hal {

class HDF5Genome;
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
   void open(const std::string& alignmentPath) const;
   void close();
   void close() const;
   void setOptionsFromParser(CLParserConstPtr parser) const;
   
   Genome* addLeafGenome(const std::string& name,
                           const std::string& parentName,
                           double branchLength);

   Genome* addRootGenome(const std::string& name,
                           double branchLength);

   void removeGenome(const std::string& name);

   Genome* insertGenome(const std::string& name,
                        const std::string& parentName,
                        const std::string& childName,
                        double upperBranchLength);

   const Genome* openGenome(const std::string& name) const;

   Genome* openGenome(const std::string& name);

   void closeGenome(const Genome* genome) const;

   std::string getRootName() const;

   std::string getParentName(const std::string& name) const;

   void updateBranchLength(const std::string& parentName,
                           const std::string& childName,
                           double length);

   double getBranchLength(const std::string& parentName,
                          const std::string& childName) const;   

   std::vector<std::string> 
   getChildNames(const std::string& name) const;

   std::vector<std::string>
   getLeafNamesBelow(const std::string& name) const;

   hal_size_t getNumGenomes() const;

   MetaData* getMetaData();

   const MetaData* getMetaData() const;

   std::string getNewickTree() const;

   std::string getVersion() const;

protected:
   // Nobody creates this class except through the interface. 
   friend AlignmentPtr hdf5AlignmentInstance();
   friend AlignmentConstPtr hdf5AlignmentInstanceReadOnly();
   friend AlignmentPtr 
   hdf5AlignmentInstance(const H5::FileCreatPropList& fileCreateProps,
                         const H5::FileAccPropList& fileAccessProps,
                         const H5::DSetCreatPropList& datasetCreateProps,
                         bool inMemory);
   friend AlignmentConstPtr 
   hdf5AlignmentInstanceReadOnly(const H5::FileCreatPropList& fileCreateProps,
                                 const H5::FileAccPropList& fileAccessProps,
                                 const H5::DSetCreatPropList& datasetCreateProps,
                                 bool inMemory);

   HDF5Alignment();
   HDF5Alignment(const H5::FileCreatPropList& fileCreateProps,
                 const H5::FileAccPropList& fileAccessProps,
                 const H5::DSetCreatPropList& datasetCreateProps,
                 bool inMemory = false);
   
   void loadTree();
   void writeTree();
   void writeVersion();
   void addGenomeToTree(const std::string& name,
                        const std::pair<std::string, double>& parentName,
                        const std::vector<std::pair<std::string, double> >&
                        childNames);


protected:

   H5::H5File* _file;
   mutable H5::FileCreatPropList _cprops;
   mutable H5::FileAccPropList _aprops;
   mutable H5::DSetCreatPropList _dcprops;
   int _flags;
   HDF5MetaData* _metaData;
   static const H5std_string MetaGroupName;
   static const H5std_string TreeGroupName;
   static const H5std_string GenomesGroupName;
   static const H5std_string VersionGroupName;
   stTree* _tree;
   mutable std::map<std::string, stTree*> _nodeMap;
   bool _dirty;
   mutable std::map<std::string, HDF5Genome*> _openGenomes;
   mutable bool _inMemory;
};

}
#endif

