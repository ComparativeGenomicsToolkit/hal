/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5TOPSEGMENT_H
#define _HDF5TOPSEGMENT_H

#include <H5Cpp.h>
#include "halTopSegment.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"

namespace hal {

class HDF5TopSegmentIterator;
class HDF5BottomSegmentIterator;

class HDF5TopSegment : public TopSegment
{
   friend class HDF5TopSegmentIterator;
   friend class HDF5BottomSegmentIterator;

public:

   /** Constructor 
    * @param genome Smart pointer to genome to which segment belongs
    * @param array HDF5 array containg segment
    * @param index Index of segment in the array */
   HDF5TopSegment(HDF5Genome* genome,
                  HDF5ExternalArray* array,
                  hal_index_t index);

   /** Destructor */
   ~HDF5TopSegment();

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

   // TOP SEGMENT INTERFACE
   hal_index_t getParentIndex() const;
   bool hasParent() const;
   void setParentIndex(hal_index_t parIdx);
   bool getParentReversed() const;
   void setParentReversed(bool isReversed);
   hal_index_t getBottomParseIndex() const;
   void setBottomParseIndex(hal_index_t botParseIdx);
   hal_offset_t getBottomParseOffset() const;
   bool hasParseDown() const;
   hal_index_t getNextParalogyIndex() const;
   bool hasNextParalogy() const;
   void setNextParalogyIndex(hal_index_t parIdx);
   hal_index_t getLeftParentIndex() const;
   hal_index_t getRightParentIndex() const;
   bool isCanonicalParalog() const;

   // HDF5 SPECIFIC
   static H5::CompType dataType();
   
private:

   static const size_t genomeIndexOffset;
   static const size_t bottomIndexOffset;
   static const size_t parIndexOffset;
   static const size_t parentIndexOffset;
   static const size_t parentReversedOffset;
   static const size_t totalSize;

   mutable HDF5ExternalArray* _array;
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};

//INLINE members
inline void HDF5TopSegment::setArrayIndex(Genome* genome, 
                                          hal_index_t arrayIndex)
{
  _genome = dynamic_cast<HDF5Genome*>(genome);
  assert(_genome != NULL);
  _array = &_genome->_topArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;
}

inline void HDF5TopSegment::setArrayIndex(const Genome* genome, 
                                          hal_index_t arrayIndex) const
{
  const HDF5Genome* h5Genome = dynamic_cast<const HDF5Genome*>(genome);
  assert(h5Genome != NULL);
  _genome = const_cast<HDF5Genome*>(h5Genome);
  _array = &_genome->_topArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;
}

inline hal_index_t HDF5TopSegment::getStartPosition() const
{
  return _array->getValue<hal_index_t>(_index, genomeIndexOffset);
}

inline hal_index_t HDF5TopSegment::getEndPosition() const
{
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t HDF5TopSegment::getLength() const
{
  return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) - 
     _array->getValue<hal_size_t>(_index, genomeIndexOffset);
}

inline const Genome* HDF5TopSegment::getGenome() const
{
  return _genome;
}

inline Genome* HDF5TopSegment::getGenome()
{
  return _genome;
}

inline const Sequence* HDF5TopSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline Sequence* HDF5TopSegment::getSequence()
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline bool HDF5TopSegment::hasParseDown() const
{
  return getBottomParseIndex() != NULL_INDEX;
}

inline hal_index_t HDF5TopSegment::getNextParalogyIndex() const
{
  return _array->getValue<hal_index_t>(_index, parIndexOffset);
}

inline bool HDF5TopSegment::hasNextParalogy() const
{
  return getNextParalogyIndex() != NULL_INDEX;
}

inline void HDF5TopSegment::setNextParalogyIndex(hal_index_t parIdx)
{
  assert(parIdx != _index);
  _array->setValue(_index, parIndexOffset, parIdx);
}

inline hal_index_t HDF5TopSegment::getParentIndex() const
{
  return _array->getValue<hal_index_t>(_index, parentIndexOffset);
}

inline bool HDF5TopSegment::hasParent() const
{
  return getParentIndex() != NULL_INDEX;
}

inline void HDF5TopSegment::setParentIndex(hal_index_t parentIndex)
{
  _array->setValue(_index, parentIndexOffset, parentIndex);
}

inline bool HDF5TopSegment::getParentReversed() const
{
  return _array->getValue<bool>(_index, parentReversedOffset); 
}

inline void HDF5TopSegment::setParentReversed(bool isReversed)
{
  _array->setValue(_index, parentReversedOffset, isReversed);
}

inline hal_index_t HDF5TopSegment::getBottomParseIndex() const
{
  return _array->getValue<hal_index_t>(_index, bottomIndexOffset);
}

inline void HDF5TopSegment::setBottomParseIndex(hal_index_t parseIndex)
{
  _array->setValue(_index, bottomIndexOffset, parseIndex);
}

inline hal_index_t HDF5TopSegment::getArrayIndex() const
{
  return _index;
}

inline bool HDF5TopSegment::leftOf(hal_index_t genomePos) const
{
  return getEndPosition() < genomePos;
}

inline bool HDF5TopSegment::rightOf(hal_index_t genomePos) const
{
  return getStartPosition() > genomePos;
}

inline bool HDF5TopSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

inline bool HDF5TopSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getTopSegmentArrayIndex();
}

inline bool HDF5TopSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == (hal_index_t)_array->getSize() - 1 || 
     _index == getSequence()->getTopSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumTopSegments() - 1;
}

inline bool HDF5TopSegment::isTop() const
{
  return true;
}

inline hal_size_t HDF5TopSegment::getMappedSegments(
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

inline hal_index_t HDF5TopSegment::getLeftParentIndex() const
{
  assert(isFirst() == false);
  HDF5TopSegment leftSeg(_genome, _array, _index - 1);
  return leftSeg.getParentIndex();
}

inline hal_index_t HDF5TopSegment::getRightParentIndex() const
{
  assert(isLast() == false);
  HDF5TopSegment rightSeg(_genome, _array, _index + 1);
  return rightSeg.getParentIndex();
}


}

#endif
