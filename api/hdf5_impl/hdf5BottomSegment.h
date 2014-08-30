/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5BOTTOMSEGMENT_H
#define _HDF5BOTTOMSEGMENT_H

#include <H5Cpp.h>
#include "halDefs.h"
#include "halBottomSegment.h"
#include "hdf5Genome.h"
#include "hdf5ExternalArray.h"

namespace hal {

class HDF5TopSegmentIterator;
class HDF5BottomSegmentIterator;

class HDF5BottomSegment : public BottomSegment
{
public:

   friend class HDF5TopSegmentIterator;
   friend class HDF5BottomSegmentIterator;

    /** Constructor 
    * @param genome Smart pointer to genome to which segment belongs
    * @param array HDF5 array containg segment
    * @param index Index of segment in the array */
   HDF5BottomSegment(HDF5Genome* genome,
                     HDF5ExternalArray* array,
                     hal_index_t index);

    /** Destructor */
   ~HDF5BottomSegment();

   // SEGMENT INTERFACE
   void setArrayIndex(Genome* genome, hal_index_t arrayIndex);
   void setArrayIndex(const Genome* genome, hal_index_t arrayIndex) const;
   const Genome* getGenome() const;
   Genome* getGenome();
   const Sequence* getSequence() const;
   Sequence* getSequence();
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;
   void getString(std::string& outString) const;
   void setCoordinates(hal_index_t startPos, hal_size_t length);
   hal_index_t getArrayIndex() const;
   bool leftOf(hal_index_t genomePos) const;
   bool rightOf(hal_index_t genomePos) const;
   bool overlaps(hal_index_t genomePos) const;
   bool isFirst() const;
   bool isLast() const;
   bool isMissingData(double nThreshold) const;
   bool isTop() const;
   hal_size_t getMappedSegments(
     std::set<MappedSegmentConstPtr>& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath,
     bool doDupes,
     hal_size_t minLength,
     const Genome *coalescenceLimit,
     const Genome *mrca) const;
   void print(std::ostream& os) const;
   
   // BOTTOM SEGMENT INTERFACE
   hal_size_t getNumChildren() const;
   hal_index_t getChildIndex(hal_size_t i) const;
   hal_index_t getChildIndexG(const Genome* childGenome) const;
   bool hasChild(hal_size_t child) const;
   bool hasChildG(const Genome* childGenome) const;
   void setChildIndex(hal_size_t i, hal_index_t childIndex);
   bool getChildReversed(hal_size_t i) const;
   void setChildReversed(hal_size_t child, bool isReversed);
   hal_index_t getTopParseIndex() const;
   void setTopParseIndex(hal_index_t parseIndex);
   hal_offset_t getTopParseOffset() const;
   bool hasParseUp() const;
   hal_index_t getLeftChildIndex(hal_size_t i) const;
   hal_index_t getRightChildIndex(hal_size_t i) const;

   // HDF5 SPECIFIC
   static H5::CompType dataType(hal_size_t numChildren);
   static hal_size_t numChildrenFromDataType(
     const H5::DataType& dataType);

private:

   static const size_t genomeIndexOffset;
   static const size_t lengthOffset;
   static const size_t topIndexOffset;
   static const size_t firstChildOffset;
   static const size_t totalSize(hal_size_t numChildren);

   mutable HDF5ExternalArray* _array;
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};


//INLINE members
inline void HDF5BottomSegment::setArrayIndex(Genome* genome, 
                                             hal_index_t arrayIndex)
{
  _genome = dynamic_cast<HDF5Genome*>(genome);
  assert(_genome != NULL);
  _array = &_genome->_bottomArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;  
}

inline void HDF5BottomSegment::setArrayIndex(const Genome* genome, 
                                             hal_index_t arrayIndex) const
{
  const HDF5Genome* h5Genome = dynamic_cast<const HDF5Genome*>(genome);
  assert(h5Genome != NULL);
  _genome = const_cast<HDF5Genome*>(h5Genome);
  _array = &_genome->_bottomArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;
}

inline hal_index_t HDF5BottomSegment::getStartPosition() const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>((hsize_t)_index, genomeIndexOffset);
}

inline hal_index_t HDF5BottomSegment::getEndPosition() const
{
  assert(_index >= 0);
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t HDF5BottomSegment::getLength() const
{
  assert(_index >= 0);
  return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) - 
     _array->getValue<hal_size_t>(_index, genomeIndexOffset);
}

inline const Genome* HDF5BottomSegment::getGenome() const
{                                               
  return _genome;
}

inline Genome* HDF5BottomSegment::getGenome()
{
  return _genome;
}

inline const Sequence* HDF5BottomSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline Sequence* HDF5BottomSegment::getSequence()
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline hal_size_t HDF5BottomSegment::getNumChildren() const
{
  return _genome->getNumChildren();
}

inline hal_index_t HDF5BottomSegment::getChildIndex(hal_size_t i) const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>(
    (hsize_t)_index, firstChildOffset + 
    i * (sizeof(hal_index_t) + sizeof(bool)));
}

inline 
hal_index_t HDF5BottomSegment::getChildIndexG(const Genome* childGenome) const
{
  assert(_index >= 0);
  return getChildIndex(_genome->getChildIndex(childGenome));
}

inline bool HDF5BottomSegment::hasChild(hal_size_t i) const
{
  return getChildIndex(i) != NULL_INDEX;
}

inline bool HDF5BottomSegment::hasChildG(const Genome* childGenome) const
{
  return getChildIndexG(childGenome) != NULL_INDEX;
}

inline void HDF5BottomSegment::setChildIndex(hal_size_t i, 
                                             hal_index_t childIndex)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, firstChildOffset + 
                   i * (sizeof(hal_index_t) + sizeof(bool)), childIndex);
}

inline bool HDF5BottomSegment::getChildReversed(hal_size_t i) const
{
  assert(_index >= 0);
  return _array->getValue<bool>(
    (hsize_t)_index, firstChildOffset + 
    i * (sizeof(hal_index_t) + sizeof(bool)) +
    sizeof(hal_index_t));
}

inline void HDF5BottomSegment::setChildReversed(hal_size_t i, 
                                                bool isReversed)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, firstChildOffset + 
                   i * (sizeof(hal_index_t) + sizeof(bool)) + 
                   sizeof(hal_index_t), isReversed);
}

inline hal_index_t HDF5BottomSegment::getTopParseIndex() const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>((hsize_t)_index, topIndexOffset);
}

inline void HDF5BottomSegment::setTopParseIndex(hal_index_t parseIndex)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, topIndexOffset, parseIndex);
}

inline bool HDF5BottomSegment::hasParseUp() const
{
  return getTopParseIndex() != NULL_INDEX;
}
  
inline hal_index_t HDF5BottomSegment::getArrayIndex() const
{
  return _index;
}

inline bool HDF5BottomSegment::leftOf(hal_index_t genomePos) const
{
  return getEndPosition() < genomePos;
}

inline bool HDF5BottomSegment::rightOf(hal_index_t genomePos) const
{
  return getStartPosition() > genomePos;
}

inline bool HDF5BottomSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

inline bool HDF5BottomSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getBottomSegmentArrayIndex();
}

inline bool HDF5BottomSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == (hal_index_t)_array->getSize() - 1 || 
     _index == getSequence()->getBottomSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumBottomSegments() - 1;
}

inline bool HDF5BottomSegment::isTop() const
{
  return false;
}

inline hal_size_t HDF5BottomSegment::getMappedSegments(
  std::set<MappedSegmentConstPtr>& outSegments,
  const Genome* tgtGenome,
  const std::set<const Genome*>* genomesOnPath,
  bool doDupes,
  hal_size_t minLength,
  const Genome *coalescenceLimit,
  const Genome *mrca) const
{
  throw hal_exception("Internal error.   HDF5 Segment interface should "
                      "at some point go through the sliced segment");
}

inline hal_index_t HDF5BottomSegment::getLeftChildIndex(hal_size_t i) const
{
  assert(isFirst() == false);
  HDF5BottomSegment leftSeg(_genome, _array, _index - 1);
  return leftSeg.getChildIndex(i);
}

inline hal_index_t HDF5BottomSegment::getRightChildIndex(hal_size_t i) const
{
  assert(isLast() == false);
  HDF5BottomSegment rightSeg(_genome, _array, _index + 1);
  return rightSeg.getChildIndex(i);
}

}


#endif
