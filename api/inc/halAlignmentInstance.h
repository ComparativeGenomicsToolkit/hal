/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTINSTANCE_H
#define _HALALIGNMENTINSTANCE_H

#include "halDefs.h"
#include "halAlignment.h"
#include "halCLParser.h"

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

/*
 * Open modes for mmap files.
 */
enum {
    MMAP_READ = 0x01,       // read-access
    MMAP_WRITE = 0x02,      // write-access
    MMAP_CREATE = 0x04,     // initialize a new file, truncate if exist
    MMAP_GROW = 0x08        // allow auto-growing, see warnings
};

/*
 * Mmap file default sizes when opening file for write access.
 */
static const size_t MMAP_DEFAULT_INIT_SIZE = 64 * GIGABYTE;
static const size_t MMAP_DEFAULT_GROW_SIZE = 64 * GIGABYTE;

/** Get an instance of an HDF5-implemented Alignment with 
 * default parameters for everything */
AlignmentPtr hdf5AlignmentInstance();
  
/** Get an instance of an HDF5-implemented Alignment from command options
 */
AlignmentPtr 
hdf5AlignmentInstance(CLParserConstPtr parser);


/** Get an instance of an HDF5-implemented Alignment while specifying 
 * @param fileCreateProps File creation properties.  fairly low-level 
 * and should generally be set to H5::FileCreatPropList::DEFAULT
 * @param fileAccessProps File access properties. Contains cache-related stuff
 * @param datasetCreateProps Compression and chunking parameters among others
 * @param inMemory Store all data in memory (overrides and disables hdf5 cache)
 */
AlignmentPtr 
hdf5AlignmentInstance(const H5::FileCreatPropList& fileCreateProps,
                      const H5::FileAccPropList& fileAccessProps,
                      const H5::DSetCreatPropList& datasetCreateProps,
                      bool inMemory = false);

/** Get read-only instance of an HDF5-implemented Alignment with
 * default parameters for everything */
AlignmentConstPtr hdf5AlignmentInstanceReadOnly();


/** Get a read-only instance of an HDF5-implemented Alignment while specifying 
 * @param fileCreateProps File creation properties.  fairly low-level 
 * and should generally be set to H5::FileCreatPropList::DEFAULT
 * @param fileAccessProps File access properties. Contains cache-related stuff
 * @param datasetCreateProps Compression and chunking parameters among others
 * @param inMemory Store all data in memory (overrides and disables hdf5 cache)
 */
AlignmentConstPtr 
hdf5AlignmentInstanceReadOnly(const H5::FileCreatPropList& fileCreateProps,
                              const H5::FileAccPropList& fileAccessProps,
                              const H5::DSetCreatPropList& datasetCreateProps,
                              bool inMemory = false);


/** Get an instance of an mmap-implemented Alignment.
 * @param fileName Path to file or URL for UDC access.
 * @param mode Access mode bit map. WARNING if MMAP_GROW is specified, pointers
 *  to file might be invalidated when allocMem() is called.
 * @param initSize Initial size to allocate when creating new (MMAP_CREATE)
 * @param growSize Addition size to allocate when growing file if MMAP_GROW is specified.
 */
AlignmentPtr 
mmapAlignmentInstance(const std::string& fileName,
                      unsigned mode = hal::MMAP_READ,
                      size_t initSize = hal::MMAP_DEFAULT_INIT_SIZE,
                      size_t growSize = hal::MMAP_DEFAULT_GROW_SIZE);

/** Get an instance of an mmap-implemented Alignment in read-only mode.
 * @param fileName Path to file or URL for UDC access.
 */
AlignmentConstPtr
mmapAlignmentInstanceReadOnly(const std::string& fileName);

    
/** Get an alignment instance from a file by automatically detecting which 
 * implementation to use.  (will currently (and probably forever more) 
 * just return an HDF5 instance since that's all that exists) 
 * @param path Path of file to open 
 * @param options Command line options information */
AlignmentPtr openHalAlignment(const std::string& path,
                              CLParserConstPtr options); 

/** Get a read-only alignment instance from a file by 
 * automatically detecting which 
 * implementation to use.  (will currently (and probably forever more) 
 * just return an HDF5 instance since that's all that exists) 
 * @param path Path of file to open 
 * @param options Command line options information */
AlignmentConstPtr openHalAlignmentReadOnly(const std::string& path,
                                           CLParserConstPtr options);

}

#endif
