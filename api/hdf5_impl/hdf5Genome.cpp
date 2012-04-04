/* Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hdf5Genome.h"

using namespace hal;
using namespace std;
using namespace H5;

HDF5Genome::HDF5Genome()
{

}


HDF5Genome::~HDF5Genome()
{

}

const std::string& HDF5Genome::getName() const
{
  return _name;
}

AlignmentPtr HDF5Genome::getAlignment()
{
  return _alignmentPtr;
}

AlignmentConstPtr HDF5Genome::getAlignment() const
{
  return _alignmentPtr;
}

hal_size_t HDF5Genome::getSequenceLength() const
{
  return _dnaArray.getSize();
}

hal_size_t HDF5Genome::getNumberTopSegments() const
{
  return _topArray.getSize();
}

hal_size_t HDF5Genome::getNumberBottomSegmentes() const
{
  return _bottomArray.getSize();
}

SegmentIteratorPtr HDF5Genome::getSegmentIterator(hal_bool_t top, 
                                               hal_index_t position)
{
  return SegmentIteratorPtr(NULL);
}

SegmentIteratorConstPtr HDF5Genome::getSegmentIterator(
     hal_bool_t top, hal_index_t position) const
{
  return SegmentIteratorConstPtr(NULL);
}
   
MetaDataPtr HDF5Genome::getMetaData()
{
  return _metaData;
}

MetaDataConstPtr HDF5Genome::getMetaData() const
{
  return _metaData;
}
  
