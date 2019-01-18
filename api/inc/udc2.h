/* udc - url data cache - a caching system that keeps blocks of data fetched from URLs in
 * sparse local files for quick use the next time the data is needed. 
 *
 * This cache is enormously simplified by there being no local _write_ to the cache,
 * just reads.  
 *
 * The overall strategy of the implementation is to have a root cache directory
 * with a subdir for each file being cached.  The directory for a single cached file
 * contains two files - "bitmap" and "sparseData" that contains information on which
 * parts of the URL are cached and the actual cached data respectively. The subdirectory name
 * associated with the file is constructed from the URL in a straightforward manner.
 *     http://genome.ucsc.edu/cgi-bin/hgGateway
 * gets mapped to:
 *     rootCacheDir/http/genome.ucsc.edu/cgi-bin/hgGateway/
 * The URL protocol is the first directory under the root, and the remainder of the
 * URL, with some necessary escaping, is used to define the rest of the cache directory
 * structure, with each '/' after the protocol line translating into another directory
 * level.
 *    
 * The bitmap file contains time stamp and size data as well as an array with one bit
 * for each block of the file that has been fetched.  Currently the block size is 8K. */

#ifndef UDC2_H
#define UDC2_H

struct udc2File;
/* Handle to a cached file.  Inside of structure mysterious unless you are udc2.c. */

struct udc2File *udc2FileMayOpen(char *url, char *cacheDir, bits32 blockSize);
/* Open up a cached file. cacheDir may be null in which case udc2DefaultDir() will be
 * used.  Return NULL if file doesn't exist.
 * If blockSize is zero, the default is used. Should be a power of two.
 */

struct udc2File *udc2FileOpen(char *url, char *cacheDir, bits32 blockSize);
/* Open up a cached file.  cacheDir may be null in which case udc2DefaultDir() will be
 * used.  Abort if if file doesn't exist.
 * If blockSize is zero, the default is used. Should be a power of two.
 */

void udc2FileClose(struct udc2File **pFile);
/* Close down cached file. */

bits64 udc2Read(struct udc2File *file, void *buf, bits64 size);
/* Read a block from file.  Return amount actually read. */

#define udc2ReadOne(file, var) udc2Read(file, &(var), sizeof(var))
/* Read one variable from file or die. */

void udc2MustRead(struct udc2File *file, void *buf, bits64 size);
/* Read a block from file.  Abort if any problem, including EOF before size is read. */

#define udc2MustReadOne(file, var) udc2MustRead(file, &(var), sizeof(var))
/* Read one variable from file or die. */

bits64 udc2ReadBits64(struct udc2File *file, boolean isSwapped);
/* Read and optionally byte-swap 64 bit entity. */

bits32 udc2ReadBits32(struct udc2File *file, boolean isSwapped);
/* Read and optionally byte-swap 32 bit entity. */

bits16 udc2ReadBits16(struct udc2File *file, boolean isSwapped);
/* Read and optionally byte-swap 16 bit entity. */

float udc2ReadFloat(struct udc2File *file, boolean isSwapped);
/* Read and optionally byte-swap floating point number. */

double udc2ReadDouble(struct udc2File *file, boolean isSwapped);
/* Read and optionally byte-swap double-precision floating point number. */

int udc2GetChar(struct udc2File *file);
/* Get next character from file or die trying. */

char *udc2ReadLine(struct udc2File *file);
/* Fetch next line from udc2 cache. */

char *udc2ReadStringAndZero(struct udc2File *file);
/* Read in zero terminated string from file.  Do a freeMem of result when done. */

char *udc2FileReadAll(char *url, char *cacheDir, size_t maxSize, size_t *retSize);
/* Read a complete file via UDC. The cacheDir may be null in which case udc2DefaultDir()
 * will be used.  If maxSize is non-zero, check size against maxSize
 * and abort if it's bigger.  Returns file data (with an extra terminal for the
 * common case where it's treated as a C string).  If retSize is non-NULL then
 * returns size of file in *retSize. Do a freeMem or freez of the returned buffer
 * when done. */

struct lineFile *udc2WrapShortLineFile(char *url, char *cacheDir, size_t maxSize);
/* Read in entire short (up to maxSize) url into memory and wrap a line file around it.
 * The cacheDir may be null in which case udc2DefaultDir() will be used.  If maxSize
 * is zero then a default value (currently 64 meg) will be used. */

void udc2Seek(struct udc2File *file, bits64 offset);
/* Seek to a particular (absolute) position in file. */

void udc2SeekCur(struct udc2File *file, bits64 offset);
/* Seek to a particular (from current) position in file. */

bits64 udc2Tell(struct udc2File *file);
/* Return current file position. */

bits64 udc2Cleanup(char *cacheDir, double maxDays, boolean testOnly);
/* Remove cached files older than maxDays old. If testOnly is set
 * no clean up is done, but the size of the files that would be
 * cleaned up is still. */

void udc2ParseUrl(char *url, char **retProtocol, char **retAfterProtocol, char **retColon);
/* Parse the URL into components that udc treats separately. 
 * *retAfterProtocol is Q-encoded to keep special chars out of filenames.  
 * Free all *ret's except *retColon when done. */

void udc2ParseUrlFull(char *url, char **retProtocol, char **retAfterProtocol, char **retColon,
		     char **retAuth);
/* Parse the URL into components that udc treats separately.
 * *retAfterProtocol is Q-encoded to keep special chars out of filenames.  
 * Free all *ret's except *retColon when done. */

char *udc2DefaultDir();
/* Get default directory for cache.  Use this for the udc2FileOpen call if you
 * don't have anything better.... */

void udc2SetDefaultDir(char *path);
/* Set default directory for cache */

void udc2DisableCache();
/* Switch off caching. Re-enable with udc2SetDefaultDir */

#define udc2DevicePrefix "udc:"
/* Prefix used by convention to indicate a file should be accessed via udc.  This is
 * followed by the local path name or a url, so in common practice you see things like:
 *     udc:http://genome.ucsc.edu/goldenPath/hg18/tracks/someTrack.bb */

struct slName *udc2FileCacheFiles(char *url, char *cacheDir);
/* Return low-level list of files used in cache. */

char *udc2PathToUrl(const char *path, char *buf, size_t size, char *cacheDir);
/* Translate path into an URL, store in buf, return pointer to buf if successful
 * and NULL if not. */

long long int udc2SizeFromCache(char *url, char *cacheDir);
/* Look up the file size from the local cache bitmap file, or -1 if there
 * is no cache for url. */

time_t udc2TimeFromCache(char *url, char *cacheDir);
/* Look up the file datetime from the local cache bitmap file, or 0 if there
 * is no cache for url. */

unsigned long udc2CacheAge(char *url, char *cacheDir);
/* Return the age in seconds of the oldest cache file.  If a cache file is
 * missing, return the current time (seconds since the epoch). */

int udc2CacheTimeout();
/* Get cache timeout (if local cache files are newer than this many seconds,
 * we won't ping the remote server to check the file size and update time). */

void udc2SetCacheTimeout(int timeout);
/* Set cache timeout (if local cache files are newer than this many seconds,
 * we won't ping the remote server to check the file size and update time). */

time_t udc2UpdateTime(struct udc2File *udc);
/* return udc->updateTime */

boolean udc2FastReadString(struct udc2File *f, char buf[256]);
/* Read a string into buffer, which must be long enough
 * to hold it.  String is in 'writeString' format. */

off_t udc2FileSize(char *url);
/* fetch remote or local file size from given URL or path */

boolean udc2Exists(char *url);
/* return true if a remote or local file exists */

boolean udc2IsLocal(char *url);
/* return true if url is not a http or ftp file, just a normal local file path */

void udc2SetLog(FILE *fp);
/* Tell UDC2 where to log its statistics. */

void udc2MMap(struct udc2File *file);
/* Enable access to underlying file as memory using mmap.  udc2MMapFetch
 * must be called to actually access regions of the file. */

void *udc2MMapFetch(struct udc2File *file, bits64 offset, bits64 size);
/* Return pointer to a region of the file in memory, ensuring that regions is
 * cached.. udc2MMap must have been called to enable access.  This must be
 * called for first access to a range of the file or erroneous (zeros) data
 * maybe returned.  Maybe called multiple times on a range or overlapping
 * returns. */


#endif /* UDC2_H */

// Local Variables:
// c-file-style: "jkent-c"
// End:
