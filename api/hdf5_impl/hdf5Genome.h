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

class Hdf5SequenceIterator;
class Hdf5Alignment;
class Hdf5Sequence;
/** 
 * HDF5 implementation of hal::Genome
 */
class Hdf5Genome : public Genome
{
   friend class Hdf5TopSegment;
   friend class Hdf5BottomSegment;
   friend class Hdf5SequenceIterator;
   friend class Hdf5Sequence;
public:

   Hdf5Genome(const std::string& name,
              Hdf5Alignment* alignment,
              H5::PortableH5Location* h5Parent,
              const H5::DSetCreatPropList& dcProps,
              bool inMemory);

   virtual ~Hdf5Genome();

    static H5::PredType dnaDataType() {
        return H5::PredType::NATIVE_UINT8;
    }

    
   // GENOME INTERFACE

   const std::string& getName() const;

   void setDimensions(
     const std::vector<hal::Sequence::Info>& sequenceDimensions,
     bool storeDNAArrays);

   void updateTopDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

   void updateBottomDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

   hal_size_t getNumSequences() const;
   
   Sequence* getSequence(const std::string& name);

   const Sequence* getSequence(const std::string& name) const;

   Sequence* getSequenceBySite(hal_size_t position);
   const Sequence* getSequenceBySite(hal_size_t position) const;
   
   SequenceIteratorPtr getSequenceIterator(
     hal_index_t position);

   SequenceIteratorPtr getSequenceIterator(
     hal_index_t position) const;

   SequenceIteratorPtr getSequenceEndIterator() const;

   MetaData* getMetaData();

   const MetaData* getMetaData() const;

   bool containsDNAArray() const;

    const Alignment* getAlignment() const;  // can't be inlined due to mutual include

    Alignment* getAlignment();  // can't be inlined due to mutual include

   void rename(const std::string &newName);

   // SEGMENTED SEQUENCE INTERFACE

   hal_size_t getSequenceLength() const;
   
   hal_size_t getNumTopSegments() const;

   hal_size_t getNumBottomSegments() const;

   TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position);

   TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position) const;

   BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position);

   BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position) const;

   DnaIteratorPtr getDnaIterator(hal_index_t position);

   DnaIteratorPtr getDnaIterator(hal_index_t position) const;

   DnaIteratorPtr getDNAEndIterator() const;

   ColumnIteratorPtr getColumnIterator(const std::set<const Genome*>* targets,
                                            hal_size_t maxInsertLength,
                                            hal_index_t position,
                                            hal_index_t lastPosition,
                                            bool noDupes,
                                            bool noAncestors,
                                            bool reverseStrand,
                                            bool unique,
                                            bool onlyOrthologs) const;

   void getString(std::string& outString) const;

   void setString(const std::string& inString);

   void getSubString(std::string& outString, hal_size_t start,
                             hal_size_t length) const;

   void setSubString(const std::string& intString, 
                             hal_size_t start,
                             hal_size_t length);

   RearrangementPtr getRearrangement(hal_index_t position,
                                     hal_size_t gapLengthThreshold,
                                     double nThreshold,
                                     bool atomic = false) const;
   
   GappedTopSegmentIteratorPtr getGappedTopSegmentIterator(
     hal_index_t i, hal_size_t gapThreshold, bool atomic) const;

   GappedBottomSegmentIteratorPtr getGappedBottomSegmentIterator(
     hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
     bool atomic) const;

   // HDF5 SPECIFIC 
   void write();
   void read();
   void create();
   void resetTreeCache();
   void resetBranchCaches();

private:

   void readSequences();
   void writeSequences(const std::vector<hal::Sequence::Info>&
                       sequenceDimensions);
   void deleteSequenceCache();
   void loadSequencePosCache() const;
   void loadSequenceNameCache() const;
   void setGenomeTopDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

   void setGenomeBottomDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);


private:

   Hdf5Alignment* _alignment;
   H5::PortableH5Location* _h5Parent;
   AlignmentPtr _alignmentPtr;
   std::string _name;
   HDF5MetaData* _metaData;
   HDF5MetaData* _rup;
   Hdf5ExternalArray _dnaArray;
   Hdf5ExternalArray _topArray;
   Hdf5ExternalArray _bottomArray;
   Hdf5ExternalArray _sequenceIdxArray;
   Hdf5ExternalArray _sequenceNameArray;
   H5::Group _group;
   H5::DSetCreatPropList _dcprops;
   hal_size_t _numChildrenInBottomArray;
   hal_size_t _totalSequenceLength;
   hal_size_t _numChunksInArrayBuffer;

   mutable std::map<hal_size_t, Hdf5Sequence*> _sequencePosCache;
   mutable std::vector<Hdf5Sequence*> _zeroLenPosCache;
   mutable std::map<std::string, Hdf5Sequence*> _sequenceNameCache;

   static const std::string dnaArrayName;
   static const std::string topArrayName;
   static const std::string bottomArrayName;
   static const std::string sequenceIdxArrayName;
   static const std::string sequenceNameArrayName;
   static const std::string metaGroupName;
   static const std::string rupGroupName;

   static const double dnaChunkScale;
};


}
#endif

