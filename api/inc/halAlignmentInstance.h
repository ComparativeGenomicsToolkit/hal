/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTINSTANCE_H
#define _HALALIGNMENTINSTANCE_H

#include "halDefs.h"
#include "halAlignment.h"


// we don't include hdf5 from our interface headers.  
namespace H5 {
class FileCreatPropList;
class FileAccPropList;
class DSetCreatPropList;
}

/** HAL API namespace */
namespace hal {

/** Get an instance of an HDF5-implemented Alignment with 
 * default parameters for everything */
AlignmentPtr hdf5AlignmentInstance();
  
/** Get read-only instance of an HDF5-implemented Alignment with
 * default parameters for everything */
AlignmentConstPtr hdf5AlignmentInstanceReadOnly();

/** Get an instance of an HDF5-implemented Alignment while specifying 
 * @param fileCreateProps File creation properties.  fairly low-level 
 * and should generally be set to H5::FileCreatPropList::DEFAULT
 * @param fileAccessProps File access properties. Contains cache-related stuff
 * @param datasetCreateProps Compression and chunking parameters among others
 */
AlignmentPtr 
hdf5AlignmentInstance(const H5::FileCreatPropList& fileCreateProps,
                      const H5::FileAccPropList& fileAccessProps,
                      const H5::DSetCreatPropList& datasetCreateProps);

/** Get a read-only instance of an HDF5-implemented Alignment while specifying 
 * @param fileCreateProps File creation properties.  fairly low-level 
 * and should generally be set to H5::FileCreatPropList::DEFAULT
 * @param fileAccessProps File access properties. Contains cache-related stuff
 * @param datasetCreateProps Compression and chunking parameters among others
 */
AlignmentConstPtr 
hdf5AlignmentInstanceReadOnly(const H5::FileCreatPropList& fileCreateProps,
                              const H5::FileAccPropList& fileAccessProps,
                              const H5::DSetCreatPropList& datasetCreateProps);

}

#endif
