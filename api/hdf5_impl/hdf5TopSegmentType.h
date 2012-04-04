/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5TOPSEGMENTTYPE_H
#define _HDF5TOPSEGMENTTYPE_H

#include "topSegmentType.h"
#includeo "rawH5ExternalArray.h"

namespace hal {

/**
 * Specialization of Top Segment type for HDF5 external array. 
 */
struct TopSegmentType<RawH5ExternalArray>
{
   typedef ArrayTraits<RawH5ExternalArray>::size_type size_type;
   typedef ArrayTraits<RawH5ExternalArray>::index_type index_type;
   typedef ArrayTraits<RawH5ExternalArray>::bool_type bool_type;
   typedef ArrayTraits<RawH5ExternalArray>::h5index_type h5index_type;
   typedef ArrayTraits<RawH5ExternalArray>::h5bool_type h5bool_type;

   /** Return object with information for arrays containing top 
    * segments */
   static H5::CompType getTypeInfo(H5::hbool_t isRoot);

   static size_type getGenomeIndexOffset();
   static size_type getLengthOffset();
   static size_type getBottomIndexOffset();
   static size_type getBottomOffsetOffset();
   static size_type getParentIndexOffset(size_type i);
   static size_type getParentInverseOffset(size_type i);
};

// INLINE MEMBERS

inline size_type BottomSegmentType<RawH5ExternalArray>::getGenomeIndexOffset()
{
  return 0;
}

inline size_type BottomSegmentType<RawH5ExternalArray>::getLengthOffset()
{
  return sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::getBottomIndexOffset()
{
  return 2 * sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::getBottomOffsetOffset()
{
  return 3 * sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::
getParentIndexOffset(size_type i)
{
  return 4 * sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::
getParentInverseOffset(size_type i)
{
  return getParentIndexOffset() + sizeof(bool_type);
}

inline H5::CompType TopSegmentType<RawH5ExternalArray>::getTypeInfo(
  hsize_t numChildren)
{
  H5::CompType dataType;
  dataType.insertMember("genomeIdx", getGenomeIndexOffset(), h5index_type);
  dataType.insertMember("length", getLengthOffset(), h5index_type);
  dataType.insertMember("bottomIdx", getBottomIndexOffset(), h5index_type);
  dataType.insertMember("bottomOffset", getBottomOffsetOffset(), h5index_type);
  if (isRoot == false)
  {
    dataType.insertMember("parentIdx", getParentIndexOffset(), h5index_type);
    dataType.insertMember("reverseFlag", getParentInverseOffset(), h5bool_type);
  }

  return dataType;
}

}
#endif
