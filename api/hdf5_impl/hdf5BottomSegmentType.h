/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5BOTTOMSEGMENTTYPE_H
#define _HDF5BOTTOMSEGMENTTYPE_H

#include <sstream>
#include "bottomSegmentType.h"
#include "rawH5ExternalArray.h"
#include "rawH5ExternalArrayTraits.h"


namespace hal {

/**
 * Specialization of Bottom Segment type for HDF5 external array. 
 */
struct BottomSegmentType<RawH5ExternalArray>
{
   typedef ArrayTraits<RawH5ExternalArray>::size_type size_type;
   typedef ArrayTraits<RawH5ExternalArray>::index_type index_type;
   typedef ArrayTraits<RawH5ExternalArray>::bool_type bool_type;
   typedef ArrayTraits<RawH5ExternalArray>::h5index_type h5index_type;
   typedef ArrayTraits<RawH5ExternalArray>::h5bool_type h5bool_type;

   /** Return object with information for arrays containing bottom
    * segments. */
   static H5::CompType getTypeInfo(size_type numChildren);
   
   static size_type getGenomeIndexOffset();
   static size_type getLengthOffset();
   static size_type getTopIndexOffset();
   static size_type getTopOffsetOffset();
   static size_type getChildIndexOffset(size_type i);
   static size_type getChildInverseOffset(size_type i);
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

inline size_type BottomSegmentType<RawH5ExternalArray>::getTopIndexOffset()
{
  return 2 * sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::getTopOffsetOffset()
{
  return 3 * sizeof(index_type);
}

inline size_type BottomSegmentType<RawH5ExternalArray>::
getChildIndexOffset(size_type i)
{
  return 3 * sizeof(index_type) + (sizeof(index_type) + sizeof(bool_type)) * i;
}

inline size_type BottomSegmentType<RawH5ExternalArray>::
getChildInverseOffset(size_type i)
{
  return getChildIndexOffset(i) + sizeof(bool_type);
}

inline H5::CompType BottomSegmentType<RawH5ExternalArray>::getTypeInfo(
  hsize_t numChildren)
{
  H5::CompType dataType;
  dataType.insertMember("genomeIdx", getGenomeIndexOffset(), h5index_type);
  dataType.insertMember("length", getLengthOffset(), h5index_type);
  dataType.insertMember("topIdx", getTopIndexOffset(), h5index_type);
  dataType.insertMember("topOffset", getTopOffsetOffset(), h5index_type);
  for(H5::hsize_t i = 0; i < numChildren; ++i)
  {
    std::stringstream ss;
    ss << i;
    H5::H5std_string number = ss.str();
    dataType.insertMember("childIdx" + number, getChildIndexOffset(i), 
                          h5index_type);
    dataType.insertMember("reverseFlag" + number, getChildInverseOffset(i),
                          h5bool_type);
  }
  return dataType;
}

}
#endif
