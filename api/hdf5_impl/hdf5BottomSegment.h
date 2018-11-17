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

class Hdf5BottomSegment : public BottomSegment
{
public:

    /** Constructor 
    * @param genome Smart pointer to genome to which segment belongs
    * @param array HDF5 array containg segment
    * @param index Index of segment in the array */
   Hdf5BottomSegment( Hdf5Genome* genome,
                     Hdf5ExternalArray* array,
                     hal_index_t index);

    /** Destructor */
   ~Hdf5BottomSegment();

   // SEGMENT INTERFACE
   void setArrayIndex(Genome* genome, hal_index_t arrayIndex);
   const Genome* getGenome() const;
   Genome* getGenome();
   const Sequence* getSequence() const;
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
   bool isTop() const;
   hal_size_t getMappedSegments(
     MappedSegmentSet& outSegments,
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

   Hdf5ExternalArray* _array;
   hal_index_t _index;
    Hdf5Genome* _genome;
};


//INLINE members
inline void Hdf5BottomSegment::setArrayIndex(Genome* genome, 
                                             hal_index_t arrayIndex)
{
  _genome = dynamic_cast<Hdf5Genome*>(genome);
  assert(_genome != NULL);
  _array = &_genome->_bottomArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;  
}

inline hal_index_t Hdf5BottomSegment::getStartPosition() const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>((hsize_t)_index, genomeIndexOffset);
}

inline hal_index_t Hdf5BottomSegment::getEndPosition() const
{
  assert(_index >= 0);
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t Hdf5BottomSegment::getLength() const
{
  assert(_index >= 0);
  return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) - 
     _array->getValue<hal_size_t>(_index, genomeIndexOffset);
}

inline const Genome* Hdf5BottomSegment::getGenome() const
{                                               
  return _genome;
}

inline Genome* Hdf5BottomSegment::getGenome()
{
    return _genome;
}

inline const Sequence* Hdf5BottomSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline hal_size_t Hdf5BottomSegment::getNumChildren() const
{
  return _genome->getNumChildren();
}

inline hal_index_t Hdf5BottomSegment::getChildIndex(hal_size_t i) const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>(
    (hsize_t)_index, firstChildOffset + 
    i * (sizeof(hal_index_t) + sizeof(bool)));
}

inline 
hal_index_t Hdf5BottomSegment::getChildIndexG(const Genome* childGenome) const
{
  assert(_index >= 0);
  return getChildIndex(_genome->getChildIndex(childGenome));
}

inline bool Hdf5BottomSegment::hasChild(hal_size_t i) const
{
  return getChildIndex(i) != NULL_INDEX;
}

inline bool Hdf5BottomSegment::hasChildG(const Genome* childGenome) const
{
  return getChildIndexG(childGenome) != NULL_INDEX;
}

inline void Hdf5BottomSegment::setChildIndex(hal_size_t i, 
                                             hal_index_t childIndex)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, firstChildOffset + 
                   i * (sizeof(hal_index_t) + sizeof(bool)), childIndex);
}

inline bool Hdf5BottomSegment::getChildReversed(hal_size_t i) const
{
  assert(_index >= 0);
  return _array->getValue<bool>(
    (hsize_t)_index, firstChildOffset + 
    i * (sizeof(hal_index_t) + sizeof(bool)) +
    sizeof(hal_index_t));
}

inline void Hdf5BottomSegment::setChildReversed(hal_size_t i, 
                                                bool isReversed)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, firstChildOffset + 
                   i * (sizeof(hal_index_t) + sizeof(bool)) + 
                   sizeof(hal_index_t), isReversed);
}

inline hal_index_t Hdf5BottomSegment::getTopParseIndex() const
{
  assert(_index >= 0);
  return _array->getValue<hal_index_t>((hsize_t)_index, topIndexOffset);
}

inline void Hdf5BottomSegment::setTopParseIndex(hal_index_t parseIndex)
{
  assert(_index >= 0);
  _array->setValue((hsize_t)_index, topIndexOffset, parseIndex);
}

inline bool Hdf5BottomSegment::hasParseUp() const
{
  return getTopParseIndex() != NULL_INDEX;
}
  
inline hal_index_t Hdf5BottomSegment::getArrayIndex() const
{
  return _index;
}

inline bool Hdf5BottomSegment::leftOf(hal_index_t genomePos) const
{
  return getEndPosition() < genomePos;
}

inline bool Hdf5BottomSegment::rightOf(hal_index_t genomePos) const
{
  return getStartPosition() > genomePos;
}

inline bool Hdf5BottomSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

inline bool Hdf5BottomSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getBottomSegmentArrayIndex();
}

inline bool Hdf5BottomSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == (hal_index_t)_array->getSize() - 1 || 
     _index == getSequence()->getBottomSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumBottomSegments() - 1;
}

inline bool Hdf5BottomSegment::isTop() const
{
  return false;
}

inline hal_size_t Hdf5BottomSegment::getMappedSegments(
  MappedSegmentSet& outSegments,
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

inline hal_index_t Hdf5BottomSegment::getLeftChildIndex(hal_size_t i) const
{
  assert(isFirst() == false);
  Hdf5BottomSegment leftSeg(_genome, _array, _index - 1);
  return leftSeg.getChildIndex(i);
}

inline hal_index_t Hdf5BottomSegment::getRightChildIndex(hal_size_t i) const
{
  assert(isLast() == false);
  Hdf5BottomSegment rightSeg(_genome, _array, _index + 1);
  return rightSeg.getChildIndex(i);
}

}


#endif
