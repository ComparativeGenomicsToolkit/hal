/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halAlignmentInstance.h"
#include "halCLParser.h"
#include "halCommon.h"
#include "hdf5Alignment.h"
#include "mmapAlignment.h"
#include <cassert>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#ifdef ENABLE_UDC
#include "udc2.h"
#endif

using namespace std;
using namespace H5;
using namespace hal;

const std::string hal::STORAGE_FORMAT_HDF5 = "hdf5";
const std::string hal::STORAGE_FORMAT_MMAP = "mmap";

/* get default FileCreatPropList with HAL default properties set */
const H5::FileCreatPropList &hal::hdf5DefaultFileCreatPropList() {
    static bool initialize = false;
    static H5::FileCreatPropList fileCreateProps;
    if (not initialize) {
        fileCreateProps.copy(H5::FileCreatPropList::DEFAULT);
        initialize = true;
    }
    return fileCreateProps;
}

/* get default FileAccPropList with HAL default properties set */
const H5::FileAccPropList &hal::hdf5DefaultFileAccPropList() {
    static bool initialize = false;
    static H5::FileAccPropList fileAccessProps;
    if (not initialize) {
#if H5_VERSION_GE(1, 14, 4)
        // hdf5 stopped working with HAL in 1.14.4
        // haven't had time to dig too deep but it seems like there's some new strictness
        // about metadata checksums or unused bits or something
        // luckily we can turn off using this (C-only) API call
        // https://support.hdfgroup.org/documentation/hdf5/latest/group___f_a_p_l.html#gafa8e677af3200e155e9208522f8e05c0        
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_relax_file_integrity_checks(fapl_id, H5F_RFIC_ALL);
        fileAccessProps.copy(H5::FileAccPropList(fapl_id));
        H5Pclose(fapl_id);
#else
        fileAccessProps.copy(H5::FileAccPropList::DEFAULT);
#endif
        fileAccessProps.setCache(Hdf5Alignment::DefaultCacheMDCElems, Hdf5Alignment::DefaultCacheRDCElems,
                                 Hdf5Alignment::DefaultCacheRDCBytes, Hdf5Alignment::DefaultCacheW0);
        initialize = true;
    }
    return fileAccessProps;
}

/* get default DSetCreatPropList  with HAL default properties set */
const H5::DSetCreatPropList &hal::hdf5DefaultDSetCreatPropList() {
    static bool initialize = false;
    static H5::DSetCreatPropList datasetCreateProps;
    if (not initialize) {
        datasetCreateProps.copy(H5::DSetCreatPropList::DEFAULT);
        datasetCreateProps.setChunk(1, &Hdf5Alignment::DefaultChunkSize);
        datasetCreateProps.setDeflate(Hdf5Alignment::DefaultCompression);
        initialize = true;
    }
    return datasetCreateProps;
}

Alignment *hal::hdf5AlignmentInstance(const std::string &alignmentPath, unsigned mode,
                                      const H5::FileCreatPropList &fileCreateProps, const H5::FileAccPropList &fileAccessProps,
                                      const H5::DSetCreatPropList &datasetCreateProps, bool inMemory) {
    return new Hdf5Alignment(alignmentPath, mode, fileCreateProps, fileAccessProps, datasetCreateProps, inMemory);
}

Alignment *hal::mmapAlignmentInstance(const std::string &alignmentPath, unsigned mode, size_t fileSize) {
    return new MMapAlignment(alignmentPath, mode, fileSize);
}

static const int DETECT_INITIAL_NUM_BYTES = 64;

static std::string udcGetInitialBytes(const std::string &path, const CLParser *options) {
#ifdef ENABLE_UDC
    struct udc2File *udcFile = udc2FileMayOpen(const_cast<char *>(path.c_str()), NULL, UDC_BLOCK_SIZE);
    if (udcFile == NULL) {
        throw hal_exception("can't open via UDC: " + path);
    }
    char buf[DETECT_INITIAL_NUM_BYTES];
    bits64 bytesRead = udc2Read(udcFile, buf, DETECT_INITIAL_NUM_BYTES);
    udc2FileClose(&udcFile);
    return string(buf, 0, bytesRead);
#else
    throw hal_exception("URL to HAL file supplied however UDC is not compiled into HAL library: " + path);
#endif
}

static std::string localGetInitialBytes(const std::string &path) {
    std::ifstream halFh;
    halFh.open(path);
    if (not halFh) {
        throw hal_errno_exception(path, "can't open HAL file", errno);
    }
    char buf[DETECT_INITIAL_NUM_BYTES];
    halFh.read(buf, DETECT_INITIAL_NUM_BYTES);
    return string(buf, 0, halFh.gcount());
}

const std::string &hal::detectHalAlignmentFormat(const std::string &path, const CLParser *options) {
    std::string initialBytes;
    if (isUrl(path)) {
        initialBytes = udcGetInitialBytes(path, options);
    } else {
        initialBytes = localGetInitialBytes(path);
    }
    if (Hdf5Alignment::isHdf5File(initialBytes)) {
        return STORAGE_FORMAT_HDF5;
    } else if (MMapFile::isMmapFile(initialBytes)) {
        return STORAGE_FORMAT_MMAP;
    } else {
        static const string empty;
        return empty;
    }
}

AlignmentPtr hal::openHalAlignment(const std::string &path, const CLParser *options, unsigned mode,
                                   const std::string &overrideFormat) {
    std::string fmt;
    if (not overrideFormat.empty()) {
        fmt = overrideFormat;
    } else if ((mode & CREATE_ACCESS) == 0) {
        fmt = detectHalAlignmentFormat(path, options);
        if (fmt.empty()) {
            throw hal_exception("unable to determine HAL storage format of " + path);
        }
    } else if (options != NULL) {
        fmt = options->getOption<const std::string &>("format");
    } else {
        fmt = STORAGE_FORMAT_HDF5;
    }
    if (fmt == STORAGE_FORMAT_HDF5) {
        if (options == NULL) {
            return AlignmentPtr(new Hdf5Alignment(path, mode, hdf5DefaultFileCreatPropList(), hdf5DefaultFileAccPropList(),
                                                  hdf5DefaultDSetCreatPropList()));
        } else {
            return AlignmentPtr(new Hdf5Alignment(path, mode, options));
        }
    } else if (fmt == STORAGE_FORMAT_MMAP) {
        if (options == NULL) {
            return AlignmentPtr(new MMapAlignment(path, mode));
        } else {
            return AlignmentPtr(new MMapAlignment(path, mode, options));
        }
    } else {
        throw hal_exception("invalid --format argument " + fmt + ", expected one of " + STORAGE_FORMAT_HDF5 + " or " +
                            STORAGE_FORMAT_MMAP);
    }
}
