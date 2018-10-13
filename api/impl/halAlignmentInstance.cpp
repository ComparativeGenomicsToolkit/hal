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
#include <sstream>
#include "halCommon.h"
#include "halAlignmentInstance.h"
#include "hdf5Alignment.h"
#include "hdf5CLParser.h"

using namespace std;
using namespace H5;
using namespace hal;

const char *hal::STORAGE_FORMAT_HDF5 = "HDF5";
const char *hal::STORAGE_FORMAT_MMAP = "mmap";

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


                           


AlignmentPtr
hal::hdf5AlignmentInstance(const std::string& alignmentPath,
                           unsigned mode,
                           const H5::FileCreatPropList& fileCreateProps,
                           const H5::FileAccPropList& fileAccessProps,
                           const H5::DSetCreatPropList& datasetCreateProps,
                           bool inMemory) {
  HDF5Alignment* al = new HDF5Alignment(alignmentPath, mode, fileCreateProps,
                                        fileAccessProps, datasetCreateProps,
                                        inMemory);
  return AlignmentPtr(al);
}

AlignmentPtr 
hal::mmapAlignmentInstance(const std::string& alignmentPath,
                           unsigned mode,
                           size_t initSize,
                           size_t growSize) {
  return AlignmentPtr(NULL);
}

AlignmentPtr hal::openHalAlignment(const std::string& path,
                                   CLParserConstPtr options,
                                   unsigned mode)
{
    /* detect which kind of file it is here (maybe by extension?) */
    
    return AlignmentPtr(new HDF5Alignment(path, mode, options));

}
