/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5GENOME_H
#define _HDF5GENOME_H

#include <H5Cpp.h>
#include "halGenome.h"
#include "hdf5ExternalArray.h"
#include "hdf5Alignment.h"
#include "halTopSegmentIterator.h"
#include "halBottomSegmentIterator.h"
#include "hdf5MetaData.h"


namespace hal {

class HDF5SegmentIterator;
class HDF5Alignment;
/** 
 * HDF5 implementation of hal::Genome
 */
class HDF5Genome : public Genome
{
   friend class HDF5SegmentIterator;
public:

   HDF5Genome(const std::string& name,
              HDF5Alignment* alignment,
              H5::CommonFG* h5Parent,
              const H5::DSetCreatPropList& dcProps);

   virtual ~HDF5Genome();

   void reset(hal_size_t totalSequenceLength,
              hal_size_t numTopSegments,
              hal_size_t numBottomSegments);

   void resetTopSegments(hal_size_t numTopSegments);

   void resetBottomSegments(hal_size_t numBottomSegments);

   const std::string& getName() const;
   AlignmentPtr getAlignment();
   hal_size_t getSequenceLength() const;
   hal_size_t getNumberTopSegments() const;
   hal_size_t getNumberBottomSegments() const;

   TopSegmentIteratorPtr getTopSegmentIterator(hal_index_t position);

   TopSegmentIteratorConstPtr getTopSegmentIterator(hal_index_t position) const;

   BottomSegmentIteratorPtr getBottomSegmentIterator(hal_index_t position);

   BottomSegmentIteratorConstPtr getBottomSegmentIterator(hal_index_t position) const;

   MetaData* getMetaData();
   const MetaData* getMetaData() const;

   Genome* getParent();
   const Genome* getParent() const;

   Genome* getChild(hal_size_t childIdx);
   const Genome* getChild(hal_size_t childIdx) const;
   hal_size_t getNumChildren() const;

   void write();
   void read();
   void create();
   void resetTreeCache();

protected:

   HDF5Alignment* _alignment;
   H5::CommonFG* _h5Parent;
   AlignmentPtr _alignmentPtr;
   std::string _name;
   HDF5MetaData* _metaData;
   HDF5ExternalArray _dnaArray;
   HDF5ExternalArray _topArray;
   HDF5ExternalArray _bottomArray;
   H5::Group _group;
   H5::DSetCreatPropList _dcprops;
   hal_size_t _numChildrenInBottomArray;

   mutable Genome* _parentCache;
   mutable std::vector<Genome*> _childCache;

   static const std::string dnaArrayName;
   static const std::string topArrayName;
   static const std::string bottomArrayName;
   static const std::string metaGroupName;
};


}
#endif

