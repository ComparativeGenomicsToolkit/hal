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

class Hdf5TopSegment : public TopSegment
{
public:

   /** Constructor 
    * @param genome Smart pointer to genome to which segment belongs
    * @param array HDF5 array containg segment
    * @param index Index of segment in the array */
   Hdf5TopSegment(Hdf5Genome* genome,
                  Hdf5ExternalArray* array,
                  hal_index_t index);

   // SEGMENT INTERFACE
   void setArrayIndex(Genome* genome, hal_index_t arrayIndex);
   const Genome* getGenome() const;
   Genome* getGenome();
   const Sequence* getSequence() const;
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;
   void setCoordinates(hal_index_t startPos, hal_size_t length);
   hal_index_t getArrayIndex() const;
   bool isFirst() const;
   bool isLast() const;
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

   Hdf5ExternalArray* _array;
   hal_index_t _index;
    Hdf5Genome* _genome;
};

//INLINE members
inline void Hdf5TopSegment::setArrayIndex(Genome* genome, 
                                          hal_index_t arrayIndex)
{
  _genome = dynamic_cast<Hdf5Genome*>(genome);
  assert(_genome != NULL);
  _array = &_genome->_topArray;
  assert(arrayIndex < (hal_index_t)_array->getSize());
  _index = arrayIndex;
}

inline hal_index_t Hdf5TopSegment::getStartPosition() const
{
  return _array->getValue<hal_index_t>(_index, genomeIndexOffset);
}

inline hal_index_t Hdf5TopSegment::getEndPosition() const
{
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t Hdf5TopSegment::getLength() const
{
  return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) - 
     _array->getValue<hal_size_t>(_index, genomeIndexOffset);
}

inline const Genome* Hdf5TopSegment::getGenome() const
{
  return _genome;
}

inline Genome* Hdf5TopSegment::getGenome()
{
  return _genome;
}

inline const Sequence* Hdf5TopSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline bool Hdf5TopSegment::hasParseDown() const
{
  return getBottomParseIndex() != NULL_INDEX;
}

inline hal_index_t Hdf5TopSegment::getNextParalogyIndex() const
{
  return _array->getValue<hal_index_t>(_index, parIndexOffset);
}

inline bool Hdf5TopSegment::hasNextParalogy() const
{
  return getNextParalogyIndex() != NULL_INDEX;
}

inline void Hdf5TopSegment::setNextParalogyIndex(hal_index_t parIdx)
{
  assert(parIdx != _index);
  _array->setValue(_index, parIndexOffset, parIdx);
}

inline hal_index_t Hdf5TopSegment::getParentIndex() const
{
  return _array->getValue<hal_index_t>(_index, parentIndexOffset);
}

inline bool Hdf5TopSegment::hasParent() const
{
  return getParentIndex() != NULL_INDEX;
}

inline void Hdf5TopSegment::setParentIndex(hal_index_t parentIndex)
{
  _array->setValue(_index, parentIndexOffset, parentIndex);
}

inline bool Hdf5TopSegment::getParentReversed() const
{
  return _array->getValue<bool>(_index, parentReversedOffset); 
}

inline void Hdf5TopSegment::setParentReversed(bool isReversed)
{
  _array->setValue(_index, parentReversedOffset, isReversed);
}

inline hal_index_t Hdf5TopSegment::getBottomParseIndex() const
{
  return _array->getValue<hal_index_t>(_index, bottomIndexOffset);
}

inline void Hdf5TopSegment::setBottomParseIndex(hal_index_t parseIndex)
{
  _array->setValue(_index, bottomIndexOffset, parseIndex);
}

inline hal_index_t Hdf5TopSegment::getArrayIndex() const
{
  return _index;
}

inline bool Hdf5TopSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getTopSegmentArrayIndex();
}

inline bool Hdf5TopSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == (hal_index_t)_array->getSize() - 1 || 
     _index == getSequence()->getTopSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumTopSegments() - 1;
}

inline hal_index_t Hdf5TopSegment::getLeftParentIndex() const
{
  assert(isFirst() == false);
  Hdf5TopSegment leftSeg(_genome, _array, _index - 1);
  return leftSeg.getParentIndex();
}

inline hal_index_t Hdf5TopSegment::getRightParentIndex() const
{
  assert(isLast() == false);
  Hdf5TopSegment rightSeg(_genome, _array, _index + 1);
  return rightSeg.getParentIndex();
}


}

#endif
// Local Variables:
// mode: c++
// End:
