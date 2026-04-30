/*
 * Copyright (C) 2024 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2024 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * LZ4 and Zstd HDF5 compression filters for HAL.
 * These are registered via H5Zregister() so they are compiled directly
 * into HAL and do not require external HDF5 plugin installation.
 */

#ifndef _HDF5FILTERS_H
#define _HDF5FILTERS_H

#include <hdf5.h>

/* Standard registered filter IDs */
#define H5Z_FILTER_LZ4  ((H5Z_filter_t)32004)
#define H5Z_FILTER_ZSTD ((H5Z_filter_t)32015)

#ifdef __cplusplus
extern "C" {
#endif

/* Register LZ4 and Zstd filters with the HDF5 library.
 * Must be called before any datasets using these filters are created or read.
 * Safe to call multiple times. */
void halRegisterCompressors(void);

#ifdef __cplusplus
}
#endif

#endif /* _HDF5FILTERS_H */
