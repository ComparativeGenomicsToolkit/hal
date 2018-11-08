/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <deque>
#include "halCommon.h"
#include "halCLParser.h"
#include "halAlignmentInstance.h"
#include "hdf5Alignment.h"
#include "mmapAlignment.h"

using namespace std;
using namespace H5;
using namespace hal;

const std::string hal::STORAGE_FORMAT_HDF5 = "hdf5";
const std::string hal::STORAGE_FORMAT_MMAP = "mmap";

/* get default FileCreatPropList with HAL default properties set */
const H5::FileCreatPropList& hal::hdf5DefaultFileCreatPropList() {
    static bool initialize = false;
    static H5::FileCreatPropList fileCreateProps;
    if (not initialize) {
        fileCreateProps.copy(H5::FileCreatPropList::DEFAULT);
        initialize = true;
    }
    return fileCreateProps;
}

/* get default FileAccPropList with HAL default properties set */
const H5::FileAccPropList& hal::hdf5DefaultFileAccPropList() {
    static bool initialize = false;
    static H5::FileAccPropList fileAccessProps;
    if (not initialize) {
        fileAccessProps.copy(H5::FileAccPropList::DEFAULT);
        fileAccessProps.setCache(HDF5Alignment::DefaultCacheMDCElems,
                                 HDF5Alignment::DefaultCacheRDCElems,
                                 HDF5Alignment::DefaultCacheRDCBytes,
                                 HDF5Alignment::DefaultCacheW0);
        initialize = true;
    }
    return fileAccessProps;
}

/* get default DSetCreatPropList  with HAL default properties set */
const H5::DSetCreatPropList& hal::hdf5DefaultDSetCreatPropList() {
    static bool initialize = false;
    static H5::DSetCreatPropList datasetCreateProps;
    if (not initialize) {
        datasetCreateProps.copy(H5::DSetCreatPropList::DEFAULT);
        datasetCreateProps.setChunk(1, &HDF5Alignment::DefaultChunkSize);
        datasetCreateProps.setDeflate(HDF5Alignment::DefaultCompression);
        initialize = true;
    }
    return datasetCreateProps;
}


                           


Alignment*
hal::hdf5AlignmentInstance(const std::string& alignmentPath,
                           unsigned mode,
                           const H5::FileCreatPropList& fileCreateProps,
                           const H5::FileAccPropList& fileAccessProps,
                           const H5::DSetCreatPropList& datasetCreateProps,
                           bool inMemory) {
  return new HDF5Alignment(alignmentPath, mode, fileCreateProps,
                           fileAccessProps, datasetCreateProps,
                           inMemory);
}

Alignment* 
hal::mmapAlignmentInstance(const std::string& alignmentPath,
                           unsigned mode,
                           size_t initSize,
                           size_t growSize) {
    return new MMapAlignment(alignmentPath, mode, initSize, growSize);
}

Alignment* hal::openHalAlignment(const std::string& path,
                                 const CLParser* options,
                                 unsigned mode,
                                 const std::string& overrideFormat)
{
    /* FIXME: detect which kind of file it is here (maybe by extension?) */
    const std::string& fmt((overrideFormat.empty()) ? options->getOption<const std::string&>("format")
                           : overrideFormat);
    if (fmt == STORAGE_FORMAT_HDF5) {
        return new HDF5Alignment(path, mode, options);
    } else if (fmt == STORAGE_FORMAT_MMAP) {
        return new MMapAlignment(path, mode, options);
    } else {
        throw hal_exception("invalid --format argument " + fmt
                            + ", expected one of " + STORAGE_FORMAT_HDF5
                            + " or " + STORAGE_FORMAT_MMAP);
    }

}
