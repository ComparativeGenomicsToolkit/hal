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

   /** Get the length of the segment (number of bases) */
   hal_size_t getLength() const;

   /** Set the length of the segment
    * @param length New length of segment */
   void setLength(hal_size_t length);

   /** Get the containing (read-only) genome */
   const Genome* getGenome() const;

   /** Get the containing genome */
   Genome* getGenome();

   /** Get the containing (read-only) sequence */
   const Sequence* getSequence() const;

   /** Get the containing sequence */
   Sequence* getSequence();

   /** Get the segment's start position in the genome */
   hal_index_t getStartPosition() const;

   /** Set the segment's start position in the genome 
    * @param startPos Start position */
   void setStartPosition(hal_index_t startPos);

   /** Get index of the next paralogous segment in the genome */
   hal_index_t getNextParalogyIndex() const;

   /** Set index of the next paralogous segment in the genome 
    * @param parIdx of next segment in same genome that is 
    * homologous to this segment */
   void setNextParalogyIndex(hal_index_t parIdx);

   /** Get flag determing if next paralogous segment aligns to the current
    * one in the reverse complement */
   hal_bool_t getNextParalogyReversed() const;

   /** Set flag determing if the next paralogous segment is reversed
    * @param parReversed flag */
   void setNextParalogyReversed(hal_bool_t parReversed);

   /** Get index of the homologous segmenet in the ancestral genome */
   hal_index_t getParentIndex() const;
   
   /** Set the index of the homologous segment in the ancestra genome 
    * @param parIdx parent index to set */
   void setParentIndex(hal_index_t parIdx);

   /** Check whether segment is mapped to parent's reverse complement */
   hal_bool_t getParentReversed() const;

   /** Set whether segment is mapped to parent's reverse complement 
    * @param isReversed Flag */
   void setParentReversed(hal_bool_t isReversed);

   /** Get the index of the bottom segment in genome that contains the
    * start coordinate of this top segment */
   hal_index_t getBottomParseIndex() const;

   /** Set the index of the bototm segment in the genome that contains the
    * start coordinate of this top segment 
    * @param botParseIndx index to set */
   void setBottomParseIndex(hal_index_t botParseIdx);

   /** Get the offset in the bottom parse segment that aligns with the
    * start coordinate of this segment */
   hal_offset_t getBottomParseOffset() const;

   /** Set the offset in the bottom parse segment that aligns with the
    * start coordinate of this segment 
    * @param botParseOffset offset */
   void setBottomParseOffset(hal_offset_t botParseOffset);

   /** Get the index of the parent of the left neighbour of this segment
    * returns NULL_INDEX if the left neighbour has no parent or the 
    * current segment is the first segment in a sequence */
   virtual hal_index_t getLeftParentIndex() const;

   /** Get the index of the parent of the right neighbour of this segment
    * returns NULL_INDEX if the right neighbour has no parent or the 
    * current segment is the first segment in a sequence */
   virtual hal_index_t getRightParentIndex() const;

   /** Test if the segment is the result of a simple inseriton (ie gap): 
    * both its left and right neighbours are adjacent in the parent
    *  (or are genome extremities) */
   bool isGapInsertion() const;
   
   /** Test if the segment is an inversion between two sets of homologous
    * segments.  ie its left and right neighbours' parents are adjacent
    * to its parent in the ancestor, but the oriernations are different */
   bool isSimpleInversion() const;

   static H5::CompType dataType();

   /** Get the index of the segment in the segment array */
   hal_index_t getArrayIndex() const;

   /** Check whether segment is the first segment of a sequence */
   bool isFirst() const;

   /** Check whether segment is the last segment of a sequence */
   bool isLast() const;

   
protected:

   static const size_t genomeIndexOffset;
   static const size_t lengthOffset;
   static const size_t bottomIndexOffset;
   static const size_t bottomOffsetOffset;
   static const size_t parIndexOffset;
   static const size_t parReversedOffset;
   static const size_t parentIndexOffset;
   static const size_t parentReversedOffset;
   static const size_t totalSize;

   mutable HDF5ExternalArray* _array;
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};

//INLINE members
inline hal_index_t HDF5TopSegment::getStartPosition() const
{
  return _array->getValue<hal_index_t>(_index, genomeIndexOffset);
}

inline void 
HDF5TopSegment::setStartPosition(hal_index_t startPos)
{
  _array->setValue(_index, genomeIndexOffset, startPos);
}

inline hal_size_t HDF5TopSegment::getLength() const
{
  return _array->getValue<hal_size_t>(_index, lengthOffset);
}

inline void 
HDF5TopSegment::setLength(hal_size_t length)
{
  _array->setValue(_index, lengthOffset, length);
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

inline hal_index_t HDF5TopSegment::getNextParalogyIndex() const
{
  return _array->getValue<hal_index_t>(_index, parIndexOffset);
}

inline void HDF5TopSegment::setNextParalogyIndex(hal_index_t parIdx)
{
  _array->setValue(_index, parIndexOffset, parIdx);
}

inline hal_bool_t HDF5TopSegment::getNextParalogyReversed() const
{
  return _array->getValue<hal_bool_t>(_index, parReversedOffset);
}

inline void HDF5TopSegment::setNextParalogyReversed(hal_bool_t parReversed)
{
  _array->setValue((hsize_t)_index, parReversedOffset, parReversed);
}

inline hal_index_t HDF5TopSegment::getParentIndex() const
{
  return _array->getValue<hal_index_t>(_index, parentIndexOffset);
}

inline void HDF5TopSegment::setParentIndex(hal_index_t parentIndex)
{
  _array->setValue(_index, parentIndexOffset, parentIndex);
}

inline hal_bool_t HDF5TopSegment::getParentReversed() const
{
  return _array->getValue<hal_bool_t>(_index, parentReversedOffset); 
}

inline void HDF5TopSegment::setParentReversed(hal_bool_t isReversed)
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
   
inline hal_offset_t HDF5TopSegment::getBottomParseOffset() const
{
  return _array->getValue<hal_offset_t>(_index, bottomOffsetOffset);
}

inline void HDF5TopSegment::setBottomParseOffset(hal_offset_t parseOffset)
{
  _array->setValue(_index, bottomOffsetOffset, parseOffset);
}

inline hal_index_t HDF5TopSegment::getArrayIndex() const
{
  return _index;
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
