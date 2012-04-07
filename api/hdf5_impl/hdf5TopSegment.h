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

class HDF5TopSegment : public TopSegment
{
public:

   /** Constructor 
    * @param genome Smart pointer to genome to which segment belongs
    * @param array HDF5 array containg segment
    * @param index Index of segment in the array */
   HDF5TopSegment(GenomePtr genome,
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
   GenomeConstPtr getGenome() const;

   /** Get the containing genome */
   GenomePtr getGenome();

   /** Get the segment's start position in the genome */
   hal_index_t getStartPosition() const;

   /** Set the segment's start position in the genome 
    * @param startPos Start position */
   void setStartPosition(hal_index_t startPos);

   /** Get a copy of the string of DNA characters associated with 
    * this segment */
   std::string getSequence() const;

   /** Get index of the next paralogous segment in the genome */
   hal_index_t getNextParalogyIndex() const;

   /** Set index of the next paralogous segment in the genome 
    * @param parIdx of next segment in same genome that is 
    * homologous to this segment */
   void setNextParalogyIndex(hal_index_t parIdx);

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

   static H5::CompType dataType();
   
protected:

   static const size_t genomeIndexOffset;
   static const size_t lengthOffset;
   static const size_t bottomIndexOffset;
   static const size_t bottomOffsetOffset;
   static const size_t parIndexOffset;
   static const size_t parentIndexOffset;
   static const size_t parentReversedOffset;
   static const size_t totalSize;

   HDF5ExternalArray* _array;
   hal_index_t _index;
   GenomePtr _genome;
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

inline GenomeConstPtr HDF5TopSegment::getGenome() const
{
  return _genome;
}

inline GenomePtr HDF5TopSegment::getGenome()
{
  return _genome;
}

inline std::string HDF5TopSegment::getSequence() const
{
  return "todo";
}

inline hal_index_t HDF5TopSegment::getNextParalogyIndex() const
{
  return _array->getValue<hal_index_t>(_index, parIndexOffset);
}

inline void HDF5TopSegment::setNextParalogyIndex(hal_index_t parIdx)
{
  _array->setValue(_index, parIndexOffset, parIdx);
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

}

#endif
