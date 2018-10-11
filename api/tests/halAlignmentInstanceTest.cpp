/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halAlignmentInstanceTest.h"
#include "halAlignmentInstance.h"
#include <H5Cpp.h>

using namespace std;
using namespace hal;

AlignmentPtr getTestAlignmentInstances(const std::string& storageFormat,
                                       const std::string& alignmentPath,
                                       unsigned mode) {
    if (storageFormat == STORAGE_FORMAT_HDF5) {
        return hdf5AlignmentInstance(alignmentPath, mode,
                                     hdf5DefaultFileCreatPropList(),
                                     hdf5DefaultFileAccPropList(),
                                     hdf5DefaultDSetCreatPropList());

    } else if (storageFormat == hal::STORAGE_FORMAT_MMAP) {
        return mmapAlignmentInstance(alignmentPath, mode);
    } else {
        throw hal_exception("invalid storage format: " + storageFormat);
    }
}
