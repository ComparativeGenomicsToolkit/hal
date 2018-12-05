/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTINSTANCE_H
#define _HALALIGNMENTINSTANCE_H

#include "halDefs.h"
#include "halAlignment.h"

/*
 * for HDF5, we don't include hdf5 from our interface headers.
 */
namespace H5 {
class FileCreatPropList;
class FileAccPropList;
class DSetCreatPropList;
}

/** HAL API namespace */
namespace hal {
    class CLParser;

/**
 * Constants defining the storage format.
 */
extern const std::string STORAGE_FORMAT_HDF5;
extern const std::string STORAGE_FORMAT_MMAP;

/*
 * Open modes for files.
 */
enum {
    READ_ACCESS = 0x01,      // read-access
    WRITE_ACCESS = 0x02,     // write-access
    CREATE_ACCESS = 0x04     // initialize a new file, truncate if exist
};

/* Default values and validate HAL mode. */
static inline unsigned halDefaultAccessMode(unsigned mode) {
    // make mode sane and validate
    if (mode & CREATE_ACCESS) {
        mode |= WRITE_ACCESS;
    }
    if (mode & WRITE_ACCESS) {
        mode |= READ_ACCESS;
    }
    if ((mode & (READ_ACCESS|WRITE_ACCESS|CREATE_ACCESS)) == 0) {
        throw hal_exception("must specify at least one of READ_ACCESS, WRITE_ACCESS, or CREATE_ACCESS on open");
    }
    return mode;
}
    
/*
 * MMap file default size when opening file for write access.
 */
static const size_t MMAP_DEFAULT_FILE_SIZE_GB = 64;
static const size_t MMAP_DEFAULT_FILE_SIZE = 64 * GIGABYTE;

/* get default FileCreatPropList with HAL default properties set */
const H5::FileCreatPropList& hdf5DefaultFileCreatPropList();

/* get default FileAccPropList with HAL default properties set */
const H5::FileAccPropList& hdf5DefaultFileAccPropList();

/* get default DSetCreatPropList  with HAL default properties set */
const H5::DSetCreatPropList& hdf5DefaultDSetCreatPropList();
    
/** Get an instance of an HDF5-implemented Alignment while specifying 
 * @param alignmentPath HDF5 file or URL
 * @param mode READ_ACCESS, WRITE_ACCESS, CREATE_ACCESS
 * @param fileCreateProps File creation properties.  Fairly low-level 
 * and should generally be set to results from hdf5DefaultFileCreatPropList().
 * @param fileAccessProps File access properties. Contains cache-related stuff.
 * Default to results from hdf5DefaultFileAccPropList().
 * @param datasetCreateProps Compression and chunking parameters among others
 * Default to results from hdf5DefaultDSetCreatPropList().
 * @param inMemory Store all data in memory (overrides and disables hdf5 cache)
 */
Alignment*
hdf5AlignmentInstance(const std::string& alignmentPath,
                      unsigned mode,
                      const H5::FileCreatPropList& fileCreateProps,
                      const H5::FileAccPropList& fileAccessProps,
                      const H5::DSetCreatPropList& datasetCreateProps,
                      bool inMemory=false);

/** Get an instance of an HDF5-implemented Alignment from command options
 */
Alignment*
hdf5AlignmentInstance(const std::string& alignmentPath,
                      unsigned mode,
                      const CLParser* parser);

/** Get an instance of an mmap-implemented Alignment.
 * @param alignmentPath Path to file or URL for UDC access.
 * @param mode Access mode bit map
 * @param fileSize Size to allocate when creating new file (CREATE_ACCESS)
 */
Alignment* 
mmapAlignmentInstance(const std::string& alignmentPath,
                      unsigned mode = hal::READ_ACCESS,
                      size_t fileSize = hal::MMAP_DEFAULT_FILE_SIZE);

/** Attempt to detect HAL alignment format, or return empty string if it doesn't 
 * appear to be a hal file */
    const std::string& detectHalAlignmentFormat(const std::string& path,
                                                const CLParser* options=NULL);

    
/** Get an alignment instance from a file by automatically detecting which 
 * implementation to use.  (will currently (and probably forever more) 
 * just return an HDF5 instance since that's all that exists) 
 * @param path Path of file to open 
 * @param options Command line options information
 * @param overrideFormat If not empty, this overrides any format in options.  Used for
 * command that both input and output HALs.
 */
Alignment* openHalAlignment(const std::string& path,
                            const CLParser* options = NULL,
                            unsigned mode = hal::READ_ACCESS,
                            const std::string& overrideFormat = ""); 
}

#endif
// Local Variables:
// mode: c++
// End:
