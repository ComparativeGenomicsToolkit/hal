/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5SEQUENCE_H
#define _HDF5SEQUENCE_H

#include <H5Cpp.h>
#include "halSequence.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"

namespace hal {

class HDF5SequenceIterator;

class HDF5Sequence : public Sequence
{
   friend class HDF5SequenceIterator;

public:

   HDF5Sequence(HDF5Genome* genome,
                HDF5ExternalArray* idxArray,
                HDF5ExternalArray* nameArray,
                hal_index_t index);

   /** Destructor */
   ~HDF5Sequence();

   // SEQUENCE INTERFACE
   std::string getName() const;

   std::string getFullName() const;

   const Genome* getGenome() const;

   Genome* getGenome();

   hal_index_t getStartPosition() const;

   hal_index_t getEndPosition() const;

   hal_index_t getArrayIndex() const;

   hal_index_t getTopSegmentArrayIndex() const;

   hal_index_t getBottomSegmentArrayIndex() const;

   // SEGMENTED SEQUENCE INTERFACE

   hal_size_t getSequenceLength() const;
   
   hal_size_t getNumTopSegments() const;

   hal_size_t getNumBottomSegments() const;

   TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position);

   TopSegmentIteratorConstPtr getTopSegmentIterator(
     hal_index_t position) const;

   TopSegmentIteratorConstPtr getTopSegmentEndIterator() const;
   
   BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position);

   BottomSegmentIteratorConstPtr getBottomSegmentIterator(
     hal_index_t position) const;

   BottomSegmentIteratorConstPtr getBottomSegmentEndIterator() const;

   DNAIteratorPtr getDNAIterator(hal_index_t position);

   DNAIteratorConstPtr getDNAIterator(hal_index_t position) const;

   DNAIteratorConstPtr getDNAEndIterator() const;

   ColumnIteratorConstPtr getColumnIterator(const std::set<const Genome*>* targets,
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
   
   GappedTopSegmentIteratorConstPtr getGappedTopSegmentIterator(
     hal_index_t i, hal_size_t gapThreshold, bool atomic) const;

   GappedBottomSegmentIteratorConstPtr getGappedBottomSegmentIterator(
     hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
     bool atomic) const;


   // LOCAL NON-INTERFACE METHODS

   static H5::CompType idxDataType();
   static H5::StrType nameDataType(hal_size_t maxNameLength);

   void set(hal_size_t startPosition, const Sequence::Info& sequenceInfo,
            hal_size_t topSegmentStartIndex,
            hal_size_t bottomSegmentStartIndex);

   void setNumTopSegments(hal_size_t numTopSegments);

   void setNumBottomSegments(hal_size_t numBottomSegments);

   void setTopSegmentArrayIndex(hal_size_t topIndex);

   void setBottomSegmentArrayIndex(hal_size_t bottomIndex); 
   
protected:

   void refreshNameCache() const;
   void refreshFullNameCache() const;
   void refreshIdxCache() const;

   static const size_t startOffset;
   static const size_t topSegmentArrayIndexOffset;
   static const size_t bottomSegmentArrayIndexOffset;
   static const size_t totalSize;
   
   mutable HDF5ExternalArray* _idxArray;
   mutable HDF5ExternalArray* _nameArray;
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};

}

#endif
