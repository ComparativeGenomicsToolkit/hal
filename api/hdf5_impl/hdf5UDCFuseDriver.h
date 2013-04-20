/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Monday, August  2, 1999
 *
 * Purpose:	The public header file for the sec2 driver.
 */
#ifndef H5FDudc_fuse_H
#define H5FDudc_fuse_H

#include "H5Ipublic.h"

#define H5FD_UDC_FUSE	(H5FD_udc_fuse_init())

extern "C" {

H5_DLL hid_t H5FD_udc_fuse_init(void);
H5_DLL void H5FD_udc_fuse_term(void);
H5_DLL herr_t H5Pset_fapl_udc_fuse(hid_t fapl_id);

void H5FD_udc_fuse_set_cache_dir(const char* cacheDir);
}

#endif
