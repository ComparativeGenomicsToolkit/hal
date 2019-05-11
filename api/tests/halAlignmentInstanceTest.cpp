/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halAlignmentInstanceTest.h"
#include "halAlignmentInstance.h"
#include <H5Cpp.h>

using namespace std;
using namespace hal;

Alignment *getTestAlignmentInstances(const std::string &storageFormat, const std::string &alignmentPath, unsigned mode) {
    if (storageFormat == STORAGE_FORMAT_HDF5) {
        return hdf5AlignmentInstance(alignmentPath, mode, hdf5DefaultFileCreatPropList(), hdf5DefaultFileAccPropList(),
                                     hdf5DefaultDSetCreatPropList());

    } else if (storageFormat == hal::STORAGE_FORMAT_MMAP) {
        // We use a default init size of only 1GiB here, because the test
        // alignments we create are relatively small.
        return mmapAlignmentInstance(alignmentPath, mode, 1024 * 1024 * 1024);
    } else {
        throw hal_exception("invalid storage format: " + storageFormat);
    }
}
