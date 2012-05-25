/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halAlignmentInstance.h"
#include "hdf5Alignment.h"

using namespace std;
using namespace H5;
using namespace hal;

AlignmentPtr hal::hdf5AlignmentInstance()
{
  return AlignmentPtr(new HDF5Alignment());
}

AlignmentConstPtr hal::hdf5AlignmentInstanceReadOnly()
{
  return AlignmentPtr(new HDF5Alignment());
}

AlignmentPtr 
hal::hdf5AlignmentInstance(const FileCreatPropList& fileCreateProps,
                           const FileAccPropList& fileAccessProps,
                           const DSetCreatPropList& datasetCreateProps)
{
  HDF5Alignment* al = new HDF5Alignment(fileCreateProps,
                                        fileAccessProps,
                                        datasetCreateProps);
  return AlignmentPtr(al);
}

AlignmentConstPtr 
hal::hdf5AlignmentInstanceReadOnly(const FileCreatPropList& fileCreateProps,
                                   const FileAccPropList& fileAccessProps,
                                   const DSetCreatPropList& datasetCreateProps)
{
  HDF5Alignment* al = new HDF5Alignment(fileCreateProps,
                                        fileAccessProps,
                                        datasetCreateProps);
  return AlignmentConstPtr(al);
}

AlignmentPtr hal::openHalAlignment(const std::string& path)
{
  // detect which kind of file it is here (maybe by extension?)
  // ...

  AlignmentPtr alignment = hdf5AlignmentInstance();
  alignment->open(path, false);
  return alignment;
}

AlignmentConstPtr hal::openHalAlignmentReadOnly(const std::string& path)
{
  // detect which kind of file it is here (maybe by extension?)
  // ...

  AlignmentConstPtr alignment = hdf5AlignmentInstanceReadOnly();
  alignment->open(path);
  return alignment;
}
