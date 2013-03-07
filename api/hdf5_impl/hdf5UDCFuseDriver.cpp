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

/* Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Wednesday, October 22, 1997
 *
 * Purpose: The C STDIO virtual file driver which only uses calls from stdio.h.
 *          This also serves as an example of coding a simple file driver,
 *          therefore, it should not use any non-public definitions.
 *
 * NOTE:    This driver is not as well tested as the standard SEC2 driver
 *          and is not intended for production use!
 */
#ifdef ENABLE_UDC

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

extern "C" {
#include "common.h"
#include "udc.h"
}
#include "hdf5.h"
#include "hdf5UDCFuseDriver.h"

#ifdef H5_HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef H5_HAVE_WIN32_API
/* The following two defines must be before any windows headers are included */
#define WIN32_LEAN_AND_MEAN    /* Exclude rarely-used stuff from Windows headers */
#define NOGDI                  /* Exclude Graphic Display Interface macros */

#include <windows.h>
#include <io.h>

/* This is not defined in the Windows header files */
#ifndef F_OK
#define F_OK 00
#endif

#endif

extern "C" {

#ifdef MAX
#undef MAX
#endif /* MAX */
#define MAX(X,Y)  ((X)>(Y)?(X):(Y))

/* The driver identification number, initialized at runtime */
static hid_t H5FD_UDC_FUSE_g = 0;

/* Flexibility to not use the default path within udc library
 * if desired */
static const char* H5FD_UDC_FUSE_CACHE_PATH = NULL;

/* The maximum number of bytes which can be written in a single I/O operation */
static size_t H5_UDC_FUSE_MAX_IO_BYTES_g = (size_t)-1;

/* File operations */
typedef enum {
    H5FD_UDC_FUSE_OP_UNKNOWN=0,
    H5FD_UDC_FUSE_OP_READ=1,
    H5FD_UDC_FUSE_OP_WRITE=2,
    H5FD_UDC_FUSE_OP_SEEK=3
} H5FD_udc_fuse_file_op;

/* The description of a file belonging to this driver. The 'eoa' and 'eof'
 * determine the amount of hdf5 address space in use and the high-water mark
 * of the file (the current size of the underlying Unix file). The 'pos'
 * value is used to eliminate file position updates when they would be a
 * no-op. Unfortunately we've found systems that use separate file position
 * indicators for reading and writing so the lseek can only be eliminated if
 * the current operation is the same as the previous operation.  When opening
 * a file the 'eof' will be set to the current file size, 'eoa' will be set
 * to zero, 'pos' will be set to H5F_ADDR_UNDEF (as it is when an error
 * occurs), and 'op' will be set to H5F_OP_UNKNOWN.
 */
typedef struct H5FD_udc_fuse_t {
    H5FD_t      pub;            /* public stuff, must be first      */
    udcFile    *ufp;            /* the file handle                  */
    int         fd;             /* file descriptor (for truncate)   */
    haddr_t     eoa;            /* end of allocated region          */
    haddr_t     eof;            /* end of file; current file size   */
    haddr_t     pos;            /* current file I/O position        */
    unsigned    write_access;   /* Flag to indicate the file was opened with write access */
    H5FD_udc_fuse_file_op op;  /* last operation */
    const char*  name;           /* path of file used for id */
#ifndef H5_HAVE_WIN32_API
    /* On most systems the combination of device and i-node number uniquely
     * identify a file.  Note that Cygwin, MinGW and other Windows POSIX
     * environments have the stat function (which fakes inodes)
     * and will use the 'device + inodes' scheme as opposed to the
     * Windows code further below.
     */
    dev_t           device;     /* file device number   */
#ifdef H5_VMS
    ino_t           inode[3];   /* file i-node number   */
#else
    ino_t           inode;      /* file i-node number   */
#endif /* H5_VMS */
#else
    /* Files in windows are uniquely identified by the volume serial
     * number and the file index (both low and high parts).
     *
     * There are caveats where these numbers can change, especially
     * on FAT file systems.  On NTFS, however, a file should keep
     * those numbers the same until renamed or deleted (though you
     * can use ReplaceFile() on NTFS to keep the numbers the same
     * while renaming).
     *
     * See the MSDN "BY_HANDLE_FILE_INFORMATION Structure" entry for
     * more information.
     *
     * http://msdn.microsoft.com/en-us/library/aa363788(v=VS.85).aspx
     */
    DWORD           nFileIndexLow;
    DWORD           nFileIndexHigh;
    DWORD           dwVolumeSerialNumber;
    
    HANDLE          hFile;      /* Native windows file handle */
#endif  /* H5_HAVE_WIN32_API */
} H5FD_udc_fuse_t;

/* Structs below specify seek / tell / truncate infterface for 
 * different platforms.  We override all with udc methods */
static int udcFseekWrapper(struct udcFile* file, long long offset, int whence)
{
  assert(whence == SEEK_SET);
  udcSeek(file, offset);
  return 0;
}

static int udcFtruncateWrapper(int, long long)
{
  assert(0);
  return -1;
}


/* Use similar structure as in H5private.h by defining Windows stuff first. */
#ifdef H5_HAVE_WIN32_API
#ifndef H5_HAVE_MINGW
    #define file_fseek      udcFseekWrapper
    #define file_offset_t   __int64
    #define file_ftruncate  udcFtruncateWrapper   /* Supported in VS 2005 or newer */
    #define file_ftell      udcTell
#endif /* H5_HAVE_MINGW */
#endif /* H5_HAVE_WIN32_API */

/* Use file_xxx to indicate these are local macros, avoiding confusing 
 * with the global HD_xxx macros. 
 * Assume fseeko, which is POSIX standard, is always supported; 
 * but prefer to use fseeko64 if supported. 
 */
#ifndef file_fseek
    #ifdef H5_HAVE_FSEEKO64
        #define file_fseek      udcFseekWrapper
        #define file_offset_t   off64_t
        #define file_ftruncate  ftruncate64
        #define file_ftell      udcTell
    #else
        #define file_fseek      udcFseekWrapper
        #define file_offset_t   off_t
        #define file_ftruncate  udcFtruncateWrapper
        #define file_ftell      udcTell
    #endif /* H5_HAVE_FSEEKO64 */
#endif /* file_fseek */

/* These macros check for overflow of various quantities.  These macros
 * assume that file_offset_t is signed and haddr_t and size_t are unsigned.
 *
 * ADDR_OVERFLOW:  Checks whether a file address of type `haddr_t'
 *      is too large to be represented by the second argument
 *      of the file seek function.
 *
 * SIZE_OVERFLOW:  Checks whether a buffer size of type `hsize_t' is too
 *      large to be represented by the `size_t' type.
 *
 * REGION_OVERFLOW:  Checks whether an address and size pair describe data
 *      which can be addressed entirely by the second
 *      argument of the file seek function.
 */
/* adding for windows NT filesystem support. */
#define MAXADDR (((haddr_t)1<<(8*sizeof(file_offset_t)-1))-1)
#define ADDR_OVERFLOW(A)  (HADDR_UNDEF==(A) || ((A) & ~(haddr_t)MAXADDR))
#define SIZE_OVERFLOW(Z)  ((Z) & ~(hsize_t)MAXADDR)
#define REGION_OVERFLOW(A,Z)  (ADDR_OVERFLOW(A) || SIZE_OVERFLOW(Z) || \
    HADDR_UNDEF==(A)+(Z) || (file_offset_t)((A)+(Z))<(file_offset_t)(A))

/* Prototypes */
static H5FD_t *H5FD_udc_fuse_open(const char *name, unsigned flags,
                 hid_t fapl_id, haddr_t maxaddr);
static herr_t H5FD_udc_fuse_close(H5FD_t *lf);
static int H5FD_udc_fuse_cmp(const H5FD_t *_f1, const H5FD_t *_f2);
static herr_t H5FD_udc_fuse_query(const H5FD_t *_f1, unsigned long *flags);
static haddr_t H5FD_udc_fuse_alloc(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, hsize_t size);
static haddr_t H5FD_udc_fuse_get_eoa(const H5FD_t *_file, H5FD_mem_t type);
static herr_t H5FD_udc_fuse_set_eoa(H5FD_t *_file, H5FD_mem_t type, haddr_t addr);
static haddr_t H5FD_udc_fuse_get_eof(const H5FD_t *_file);
static herr_t  H5FD_udc_fuse_get_handle(H5FD_t *_file, hid_t fapl, void** file_handle);
static herr_t H5FD_udc_fuse_read(H5FD_t *lf, H5FD_mem_t type, hid_t fapl_id, haddr_t addr,
                size_t size, void *buf);
static herr_t H5FD_udc_fuse_write(H5FD_t *lf, H5FD_mem_t type, hid_t fapl_id, haddr_t addr,
                size_t size, const void *buf);
static herr_t H5FD_udc_fuse_flush(H5FD_t *_file, hid_t dxpl_id, unsigned closing);
static herr_t H5FD_udc_fuse_truncate(H5FD_t *_file, hid_t dxpl_id, hbool_t closing);

static const H5FD_class_t H5FD_udc_fuse_g = {
    "udc_fuse",                    /* name         */
    MAXADDR,                    /* maxaddr      */
    H5F_CLOSE_WEAK,             /* fc_degree    */
    NULL,                       /* sb_size      */
    NULL,                       /* sb_encode    */
    NULL,                       /* sb_decode    */
    0,                          /* fapl_size    */
    NULL,                       /* fapl_get     */
    NULL,                       /* fapl_copy    */
    NULL,                       /* fapl_free    */
    0,                          /* dxpl_size    */
    NULL,                       /* dxpl_copy    */
    NULL,                       /* dxpl_free    */
    H5FD_udc_fuse_open,            /* open         */
    H5FD_udc_fuse_close,           /* close        */
    H5FD_udc_fuse_cmp,             /* cmp          */
    H5FD_udc_fuse_query,           /* query        */
    NULL,                       /* get_type_map */
    H5FD_udc_fuse_alloc,           /* alloc        */
    NULL,                       /* free         */
    H5FD_udc_fuse_get_eoa,         /* get_eoa      */
    H5FD_udc_fuse_set_eoa,         /* set_eoa      */
    H5FD_udc_fuse_get_eof,         /* get_eof      */
    H5FD_udc_fuse_get_handle,      /* get_handle   */
    H5FD_udc_fuse_read,            /* read         */
    H5FD_udc_fuse_write,           /* write        */
    H5FD_udc_fuse_flush,           /* flush        */
    H5FD_udc_fuse_truncate,        /* truncate     */
    NULL,                       /* lock         */
    NULL,                       /* unlock       */
    H5FD_FLMAP_SINGLE           /* fl_map       */
};

void H5FD_udc_fuse_set_cache_dir(const char* cacheDir)
{
  H5FD_UDC_FUSE_CACHE_PATH = cacheDir;
}
/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_init
 *
 * Purpose:  Initialize this driver by registering the driver with the
 *    library.
 *
 * Return:  Success:  The driver ID for the udc_fuse driver.
 *
 *    Failure:  Negative.
 *
 * Programmer:  Robb Matzke
 *              Thursday, July 29, 1999
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5FD_udc_fuse_init(void)
{
    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    if (H5I_VFL!=H5Iget_type(H5FD_UDC_FUSE_g))
        H5FD_UDC_FUSE_g = H5FDregister(&H5FD_udc_fuse_g);

    return H5FD_UDC_FUSE_g;
} /* end H5FD_udc_fuse_init() */


/*---------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_term
 *
 * Purpose:  Shut down the VFD
 *
 * Returns:     None
 *
 * Programmer:  Quincey Koziol
 *              Friday, Jan 30, 2004
 *
 *---------------------------------------------------------------------------
 */
void
H5FD_udc_fuse_term(void)
{
    /* Reset VFL ID */
    H5FD_UDC_FUSE_g = 0;

    return;
} /* end H5FD_udc_fuse_term() */


/*-------------------------------------------------------------------------
 * Function:  H5Pset_fapl_udc_fuse
 *
 * Purpose:  Modify the file access property list to use the H5FD_UDC_FUSE
 *    driver defined in this source file.  There are no driver
 *    specific properties.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *    Thursday, February 19, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_fapl_udc_fuse(hid_t fapl_id)
{
    static const char *func = "H5FDset_fapl_udc_fuse";  /*for error reporting*/

    /*NO TRACE*/

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    if(0 == H5Pisa_class(fapl_id, H5P_FILE_ACCESS))
        H5Epush_ret(func, H5E_ERR_CLS, H5E_PLIST, H5E_BADTYPE, "not a file access property list", -1)

    return H5Pset_driver(fapl_id, H5FD_UDC_FUSE, NULL);
} /* end H5Pset_fapl_udc_fuse() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_open
 *
 * Purpose:  Create and/or opens a Standard C file as an HDF5 file.
 *
 * Errors:
 *  IO  CANTOPENFILE    File doesn't exist and CREAT wasn't
 *                      specified.
 *  IO  CANTOPENFILE    fopen() failed.
 *  IO  FILEEXISTS      File exists but CREAT and EXCL were
 *                      specified.
 *
 * Return:
 *      Success:    A pointer to a new file data structure. The
 *                  public fields will be initialized by the
 *                  caller, which is always H5FD_open().
 *
 *      Failure:    NULL
 *
 * Programmer:  Robb Matzke
 *    Wednesday, October 22, 1997
 *
 *-------------------------------------------------------------------------
 */
static H5FD_t *
H5FD_udc_fuse_open( const char *name, unsigned flags, hid_t fapl_id,
    haddr_t maxaddr)
{
    udcFile                *f = NULL;
    unsigned            write_access = 0;           /* File opened with write access? */
    H5FD_udc_fuse_t        *file = NULL;
    static const char   *func = "H5FD_udc_fuse_open";  /* Function Name for error reporting */
#ifdef H5_HAVE_WIN32_API
    struct _BY_HANDLE_FILE_INFORMATION fileinfo;
#else /* H5_HAVE_WIN32_API */
    struct stat         sb;
#endif  /* H5_HAVE_WIN32_API */
    /* Sanity check on file offsets */
    assert(sizeof(file_offset_t) >= sizeof(size_t));

    /* Quiet compiler */
    fapl_id = fapl_id;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Check arguments */
    if (!name || !*name)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_ARGS, H5E_BADVALUE, "invalid file name", NULL)
    if (0 == maxaddr || HADDR_UNDEF == maxaddr)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_ARGS, H5E_BADRANGE, "bogus maxaddr", NULL)
    if (ADDR_OVERFLOW(maxaddr))
        H5Epush_ret(func, H5E_ERR_CLS, H5E_ARGS, H5E_OVERFLOW, "maxaddr too large", NULL)

    /* Attempt to open/create the file */
        
      f = udcFileMayOpen((char*)name, (char*)H5FD_UDC_FUSE_CACHE_PATH);
    
    if (!f)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_CANTOPENFILE, "fopen failed", NULL)

    /* Build the return value */
    if(NULL == (file = (H5FD_udc_fuse_t *)calloc((size_t)1, sizeof(H5FD_udc_fuse_t)))) {
        udcFileClose(&f);
        H5Epush_ret(func, H5E_ERR_CLS, H5E_RESOURCE, H5E_NOSPACE, "memory allocation failed", NULL)
    } /* end if */
    file->ufp = f;
    file->name = name;
    file->op = H5FD_UDC_FUSE_OP_SEEK;
    file->pos = HADDR_UNDEF;
    file->write_access = write_access;    /* Note the write_access for later */
    /* note -- do we add interface to modify cache dir? */

    long long int udcSizeVal = udcSizeFromCache((char*)name, 
                                                (char*)H5FD_UDC_FUSE_CACHE_PATH);
    file->eof = udcSizeVal;
    /* everything about udc cache works for files and urls but the above, 
     * which only works for ursl.  if it fails, we try as a file*/
    if (udcSizeVal < 0)
    {
      FILE* tempHandle = fopen(name, "r");
      if (tempHandle)
      {
        fseek(tempHandle, 0, SEEK_END);
        file->eof = (haddr_t) ftell(tempHandle);
        fclose(tempHandle);
      }
    }

    /* Get the file descriptor (needed for truncate and some Windows information) */
    file->fd = 0;
    if(file->fd < 0)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_FILE, H5E_CANTOPENFILE, "unable to get file descriptor", NULL);



    return (H5FD_t*)file;
} /* end H5FD_udc_fuse_open() */


/*-------------------------------------------------------------------------
 * Function:  H5F_udc_fuse_close
 *
 * Purpose:  Closes a file.
 *
 * Errors:
 *    IO    CLOSEERROR  Fclose failed.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *    Wednesday, October 22, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_close(H5FD_t *_file)
{
    H5FD_udc_fuse_t  *file = (H5FD_udc_fuse_t*)_file;
    static const char *func = "H5FD_udc_fuse_close";  /* Function Name for error reporting */

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);
    
    udcFileClose(&file->ufp);
    
    return 0;
} /* end H5FD_udc_fuse_close() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_cmp
 *
 * Purpose:  Compares two files belonging to this driver using an
 *    arbitrary (but consistent) ordering.
 *
 * Return:
 *      Success:    A value like strcmp()
 *
 *      Failure:    never fails (arguments were checked by the caller).
 *
 * Programmer:  Robb Matzke
 *              Thursday, July 29, 1999
 *
 *-------------------------------------------------------------------------
 */
static int
H5FD_udc_fuse_cmp(const H5FD_t *_f1, const H5FD_t *_f2)
{
    const H5FD_udc_fuse_t  *f1 = (const H5FD_udc_fuse_t*)_f1;
    const H5FD_udc_fuse_t  *f2 = (const H5FD_udc_fuse_t*)_f2;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    assert(f1->name and f2->name);
    return strcmp(f1->name, f2->name);
} /* H5FD_udc_fuse_cmp() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_query
 *
 * Purpose:  Set the flags that this VFL driver is capable of supporting.
 *              (listed in H5FDpublic.h)
 *
 * Return:  Success:  non-negative
 *
 *    Failure:  negative
 *
 * Programmer:  Quincey Koziol
 *              Friday, August 25, 2000
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_query(const H5FD_t *_f, unsigned long *flags /* out */)
{
    /* Quiet the compiler */
    _f=_f;

    /* Set the VFL feature flags that this driver supports */
    if(flags) {
        *flags = 0;
        *flags|=H5FD_FEAT_AGGREGATE_METADATA; /* OK to aggregate metadata allocations */
        *flags|=H5FD_FEAT_ACCUMULATE_METADATA; /* OK to accumulate metadata for faster writes */
        *flags|=H5FD_FEAT_DATA_SIEVE;       /* OK to perform data sieving for faster raw data reads & writes */
        *flags|=H5FD_FEAT_AGGREGATE_SMALLDATA; /* OK to aggregate "small" raw data allocations */
    }

    return 0;
} /* end H5FD_udc_fuse_query() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_alloc
 *
 * Purpose:     Allocates file memory. If fseeko isn't available, makes
 *              sure the file size isn't bigger than 2GB because the
 *              parameter OFFSET of fseek is of the type LONG INT, limiting
 *              the file size to 2GB.
 *
 * Return:
 *      Success:    Address of new memory
 *
 *      Failure:    HADDR_UNDEF
 *
 * Programmer:  Raymond Lu
 *              30 March 2007
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD_udc_fuse_alloc(H5FD_t *_file, H5FD_mem_t /*UNUSED*/ type, hid_t /*UNUSED*/ dxpl_id, hsize_t size)
{
    H5FD_udc_fuse_t    *file = (H5FD_udc_fuse_t*)_file;
    haddr_t         addr;

    /* Quiet compiler */
    type = type;
    dxpl_id = dxpl_id;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Compute the address for the block to allocate */
    addr = file->eoa;

    /* Check if we need to align this block */
    if(size >= file->pub.threshold) {
        /* Check for an already aligned block */
        if((addr % file->pub.alignment) != 0)
            addr = ((addr / file->pub.alignment) + 1) * file->pub.alignment;
    } /* end if */

    file->eoa = addr + size;

    return addr;
} /* end H5FD_udc_fuse_alloc() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_get_eoa
 *
 * Purpose:  Gets the end-of-address marker for the file. The EOA marker
 *           is the first address past the last byte allocated in the
 *           format address space.
 *
 * Return:  Success:  The end-of-address marker.
 *
 *    Failure:  HADDR_UNDEF
 *
 * Programmer:  Robb Matzke
 *              Monday, August  2, 1999
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD_udc_fuse_get_eoa(const H5FD_t *_file, H5FD_mem_t /*UNUSED*/ type)
{
    const H5FD_udc_fuse_t *file = (const H5FD_udc_fuse_t *)_file;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Quiet compiler */
    type = type;

    return file->eoa;
} /* end H5FD_udc_fuse_get_eoa() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_set_eoa
 *
 * Purpose:  Set the end-of-address marker for the file. This function is
 *    called shortly after an existing HDF5 file is opened in order
 *    to tell the driver where the end of the HDF5 data is located.
 *
 * Return:  Success:  0
 *
 *    Failure:  Does not fail
 *
 * Programmer:  Robb Matzke
 *              Thursday, July 29, 1999
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_set_eoa(H5FD_t *_file, H5FD_mem_t /*UNUSED*/ type, haddr_t addr)
{
    H5FD_udc_fuse_t  *file = (H5FD_udc_fuse_t*)_file;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Quiet the compiler */
    type = type;

    file->eoa = addr;

    return 0;
}


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_get_eof
 *
 * Purpose:  Returns the end-of-file marker, which is the greater of
 *    either the Unix end-of-file or the HDF5 end-of-address
 *    markers.
 *
 * Return:  Success:  End of file address, the first address past
 *        the end of the "file", either the Unix file
 *        or the HDF5 file.
 *
 *    Failure:  HADDR_UNDEF
 *
 * Programmer:  Robb Matzke
 *              Thursday, July 29, 1999
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD_udc_fuse_get_eof(const H5FD_t *_file)
{
    const H5FD_udc_fuse_t  *file = (const H5FD_udc_fuse_t *)_file;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    return MAX(file->eof, file->eoa);
} /* end H5FD_udc_fuse_get_eof() */


/*-------------------------------------------------------------------------
 * Function:       H5FD_udc_fuse_get_handle
 *
 * Purpose:        Returns the file handle of udc_fuse file driver.
 *
 * Returns:        Non-negative if succeed or negative if fails.
 *
 * Programmer:     Raymond Lu
 *                 Sept. 16, 2002
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_get_handle(H5FD_t *_file, hid_t fapl, void** file_handle)
{
    H5FD_udc_fuse_t       *file = (H5FD_udc_fuse_t *)_file;
    static const char  *func = "H5FD_udc_fuse_get_handle";  /* Function Name for error reporting */

    /* Quiet the compiler */
    fapl = fapl;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    *file_handle = &(file->ufp);
    if(*file_handle == NULL)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_WRITEERROR, "get handle failed", -1)

    return 0;
} /* end H5FD_udc_fuse_get_handle() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_read
 *
 * Purpose:  Reads SIZE bytes beginning at address ADDR in file LF and
 *    places them in buffer BUF.  Reading past the logical or
 *    physical end of file returns zeros instead of failing.
 *
 * Errors:
 *    IO    READERROR  fread failed.
 *    IO    SEEKERROR  fseek failed.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *    Wednesday, October 22, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_read(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr, size_t size,
    void *buf/*out*/)
{
    H5FD_udc_fuse_t    *file = (H5FD_udc_fuse_t*)_file;
    static const char *func = "H5FD_udc_fuse_read";  /* Function Name for error reporting */

    /* Quiet the compiler */
    type = type;
    dxpl_id = dxpl_id;

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    /* Check for overflow */
    if (HADDR_UNDEF==addr)
        H5Epush_ret (func, H5E_ERR_CLS, H5E_IO, H5E_OVERFLOW, "file address overflowed", -1)
    if (REGION_OVERFLOW(addr, size))
        H5Epush_ret (func, H5E_ERR_CLS, H5E_IO, H5E_OVERFLOW, "file address overflowed", -1)
    if((addr + size) > file->eoa)
        H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_OVERFLOW, "file address overflowed", -1)

    /* Check easy cases */
    if (0 == size)
        return 0;
    if ((haddr_t)addr >= file->eof) {
        memset(buf, 0, size);
        return 0;
    }

    /* Seek to the correct file position. */
    if (!(file->op == H5FD_UDC_FUSE_OP_READ || file->op == H5FD_UDC_FUSE_OP_SEEK) ||
            file->pos != addr) {
        if (file_fseek(file->ufp, (file_offset_t)addr, SEEK_SET) < 0) {
            file->op = H5FD_UDC_FUSE_OP_UNKNOWN;
            file->pos = HADDR_UNDEF;
            H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_SEEKERROR, "fseek failed", -1)
        }
        file->pos = addr;
    }

    /* Read zeros past the logical end of file (physical is handled below) */
    if (addr + size > file->eof) {
        size_t nbytes = (size_t) (addr + size - file->eof);
        memset((unsigned char *)buf + size - nbytes, 0, nbytes);
        size -= nbytes;
    }

    /* Read the data.  Since we're reading single-byte values, a partial read
     * will advance the file position by N.  If N is zero or an error
     * occurs then the file position is undefined.
     */
    while(size > 0) {

        size_t bytes_in        = 0;    /* # of bytes to read       */
        size_t bytes_read      = 0;    /* # of bytes actually read */
        size_t item_size       = 1;    /* size of items in bytes */

        if(size > H5_UDC_FUSE_MAX_IO_BYTES_g)
            bytes_in = H5_UDC_FUSE_MAX_IO_BYTES_g;
        else
            bytes_in = size;

        bytes_read = udcRead(file->ufp, buf, item_size * bytes_in);

        /*
        if(0 == bytes_read && ferror(file->ufp)) { 
            file->op = H5FD_UDC_FUSE_OP_UNKNOWN;
            file->pos = HADDR_UNDEF;
            H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_READERROR, "fread failed", -1)
        }*/
        
        if(0 == bytes_read && udcTell(file->ufp) >= file->eof) {
            /* end of file but not end of format address space */
            memset((unsigned char *)buf, 0, size);
            break;
        } /* end if */
        
        size -= bytes_read;
        addr += (haddr_t)bytes_read;
        buf = (char *)buf + bytes_read;
    } /* end while */

    /* Update the file position data. */
    file->op = H5FD_UDC_FUSE_OP_READ;
    file->pos = addr;

    return 0;
}


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_write
 *
 * Purpose:  Writes SIZE bytes from the beginning of BUF into file LF at
 *    file address ADDR.
 *
 * Errors:
 *    IO    SEEKERROR   fseek failed.
 *    IO    WRITEERROR  fwrite failed.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *    Wednesday, October 22, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_write(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr,
    size_t size, const void *buf)
{
    H5FD_udc_fuse_t    *file = (H5FD_udc_fuse_t*)_file;
    static const char *func = "H5FD_udc_fuse_write";  /* Function Name for error reporting */

    /* Clear the error stack */
    H5Eclear2(H5E_DEFAULT);

    H5Epush_ret (func, H5E_ERR_CLS, H5E_IO, H5E_WRITEERROR, "udc driver cannot write!", -1)
    
    return -1;
}


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_flush
 *
 * Purpose:  Makes sure that all data is on disk.
 *
 * Errors:
 *    IO    SEEKERROR     fseek failed.
 *    IO    WRITEERROR    fflush or fwrite failed.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *    Wednesday, October 22, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_flush(H5FD_t *_file, hid_t dxpl_id, unsigned closing)
{
    H5FD_udc_fuse_t  *file = (H5FD_udc_fuse_t*)_file;
    static const char *func = "H5FD_udc_fuse_flush";  /* Function Name for error reporting */
    H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_WRITEERROR, "udc driver cannot flush", -1)
    
    return -1;
} /* end H5FD_udc_fuse_flush() */


/*-------------------------------------------------------------------------
 * Function:  H5FD_udc_fuse_truncate
 *
 * Purpose:  Makes sure that the true file size is the same (or larger)
 *    than the end-of-address.
 *
 * Errors:
 *    IO    SEEKERROR     fseek failed.
 *    IO    WRITEERROR    fflush or fwrite failed.
 *
 * Return:  Non-negative on success/Negative on failure
 *
 * Programmer:  Quincey Koziol
 *    Thursday, January 31, 2008
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_udc_fuse_truncate(H5FD_t *_file, hid_t dxpl_id, hbool_t closing)
{
    H5FD_udc_fuse_t  *file = (H5FD_udc_fuse_t*)_file;
    static const char *func = "H5FD_udc_fuse_truncate";  /* Function Name for error reporting */
    
    H5Epush_ret(func, H5E_ERR_CLS, H5E_IO, H5E_WRITEERROR, "udc driver cannot truncate", -1)
    
    return 0;
} /* end H5FD_udc_fuse_truncate() */


#ifdef _H5private_H
/*
 * This is not related to the functionality of the driver code.
 * It is added here to trigger warning if HDF5 private definitions are included
 * by mistake.  The code should use only HDF5 public API and definitions.
 */
#error "Do not use HDF5 private definitions"
#endif

}

#endif

