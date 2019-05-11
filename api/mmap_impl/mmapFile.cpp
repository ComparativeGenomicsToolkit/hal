#include "mmapFile.h"
#include "halCommon.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef ENABLE_UDC
extern "C" {
#include "common.h"
#include "udc2.h"
}
#endif

/* constants for header */
static const std::string FORMAT_NAME = "HAL-MMAP";
static const std::string MMAP_VERSION = "1.0";

/* check if first bit of file has MMAP header */
bool hal::MMapFile::isMmapFile(const std::string &initialBytes) {
    return initialBytes.compare(0, FORMAT_NAME.size(), FORMAT_NAME) == 0;
}

/* constructor, used only by derived classes */
hal::MMapFile::MMapFile(const std::string alignmentPath, unsigned mode, bool mustFetch)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _basePtr(NULL), _fileSize(0), _mustFetch(mustFetch) {
}

/* error if file is not open for write accecss */
void hal::MMapFile::validateWriteAccess() const {
    if ((_mode & WRITE_ACCESS) == 0) {
        throw hal_exception(_alignmentPath + " is not open for write access");
    }
}

/* setup pointer to header */
void hal::MMapFile::setHeaderPtr() {
    fetchIfNeeded(0, sizeof(MMapHeader));
    _header = static_cast<MMapHeader *>(_basePtr);
}

/* validate the file header and save a pointer to it. */
void hal::MMapFile::loadHeader(bool markDirty) {
    if (_fileSize < sizeof(MMapHeader)) {
        throw hal_exception(_alignmentPath + ": file size of " + std::to_string(_fileSize) + " is less that header size of " +
                            std::to_string(sizeof(MMapHeader)));
    }
    setHeaderPtr();

    // don't print found strings as it might have garbage
    if (::strcmp(_header->format, FORMAT_NAME.c_str()) != 0) {
        throw hal_exception(_alignmentPath + ": invalid file header, expected format name of '" + FORMAT_NAME + "'");
    }
    // FIXME: to need to check version compatibility
    if (::strcmp(_header->mmapVersion, MMAP_VERSION.c_str()) != 0) {
        throw hal_exception(_alignmentPath + ": incompatible mmap format versions: " + "file version " + _header->mmapVersion +
                            ", mmap API version " + MMAP_VERSION);
    }
    if ((_header->nextOffset < sizeof(MMapHeader)) or (_header->nextOffset > _fileSize)) {
        throw hal_exception(_alignmentPath + ": header nextOffset field out of bounds, probably file corruption");
    }
    if ((_header->rootOffset < sizeof(MMapHeader)) or (_header->rootOffset > _fileSize)) {
        throw hal_exception(_alignmentPath + ": header rootOffset field out of bounds, probably file corruption");
    }
    if (_header->dirty) {
        throw hal_exception(_alignmentPath + ": file is marked as dirty, most likely an inconsistent state.");
    }
    if (markDirty) {
        _header->dirty = true;
    }
}

/* create the header */
void hal::MMapFile::createHeader() {
    assert(_mode & WRITE_ACCESS);
    setHeaderPtr();
    assert(FORMAT_NAME.size() < sizeof(_header->format));
    strncpy(_header->format, FORMAT_NAME.c_str(), sizeof(_header->format) - 1);
    assert(MMAP_VERSION.size() < sizeof(_header->mmapVersion));
    strncpy(_header->mmapVersion, MMAP_VERSION.c_str(), sizeof(_header->mmapVersion) - 1);
    assert(HAL_VERSION.size() < sizeof(_header->halVersion));
    strncpy(_header->halVersion, HAL_VERSION.c_str(), sizeof(_header->halVersion) - 1);
    _header->nextOffset = alignRound(sizeof(MMapHeader));
    _header->dirty = true;
    _header->nextOffset = _header->nextOffset;
}

namespace hal {
    /* Class that implements local file version of MMapFile */
    class MMapFileLocal : public MMapFile {
      public:
        MMapFileLocal(const std::string &alignmentPath, unsigned mode, size_t fileSize);
        virtual void close();
        virtual ~MMapFileLocal();
        virtual bool isUdcProtocol() const {
            return false;
        }

      private:
        int openFile();
        void closeFile();
        void adjustFileSize(size_t size);
        void *mapFile(void *requiredAddr = NULL);
        void unmapFile();
        void openRead();
        void openWrite(size_t fileSize);

        int _fd; // open file descriptor
    };
}

/* Constructor. Open or create the specified file. */
hal::MMapFileLocal::MMapFileLocal(const std::string &alignmentPath, unsigned mode, size_t fileSize)
    : MMapFile(alignmentPath, mode, false), _fd(-1) {
    if (_mode & WRITE_ACCESS) {
        openWrite(fileSize);
    } else {
        openRead();
    }
}

/* close file, marking as clean.  Don't  */
void hal::MMapFileLocal::close() {
    if (_basePtr == NULL) {
        throw hal_exception(_alignmentPath + ": MMapFile::close() called on closed file");
    }
    if (_mode & WRITE_ACCESS) {
        adjustFileSize(_header->nextOffset);
        _header->dirty = false;
    }
    unmapFile();
    closeFile();
}

/* Destructor. write fields to header and close.  If write access and close
 * has not been called, file will me left mark dirty */
hal::MMapFileLocal::~MMapFileLocal() {
    unmapFile();
    closeFile();
}

/* open the file for the specified mode */
int hal::MMapFileLocal::openFile() {
    assert(_fd < 0);
    unsigned openMode = 0;
    if (_mode & WRITE_ACCESS) {
        openMode = O_RDWR | ((_mode & CREATE_ACCESS) ? (O_CREAT | O_TRUNC) : 0);
    } else {
        openMode = O_RDONLY;
    }
    int fd = ::open(_alignmentPath.c_str(), openMode, 0666);
    if (fd < 0) {
            throw hal_errno_exception(_alignmentPath, "open failed", errno);
    }
    return fd;
}

/* change size size of the file, possibly deleting data. */
void hal::MMapFileLocal::adjustFileSize(size_t size) {
try_truncate:
    if (ftruncate(_fd, size) < 0) {
        throw hal_errno_exception(_alignmentPath, "set size failed", errno);
    }
    _fileSize = size;
}

/* map file into memory */
void *hal::MMapFileLocal::mapFile(void *requiredAddr) {
    assert(_basePtr == NULL);
    unsigned prot = PROT_READ | ((_mode & WRITE_ACCESS) ? PROT_WRITE : 0);
    int flags = MAP_SHARED | MAP_FILE;
    if (requiredAddr != NULL) {
        // We don't want MAP_FIXED when we don't have an address we
        // need, as that will, apparently, happily map NULL to the
        // beginning of our file. That could mask some bugs. But we
        // *do* want MAP_FIXED otherwise.
        flags |= MAP_FIXED;
    }
    void *ptr = mmap(requiredAddr, _fileSize, prot, flags, _fd, 0);
    if (ptr == MAP_FAILED) {
        throw hal_errno_exception(_alignmentPath, "mmap failed", errno);
    }
    if (requiredAddr != NULL && ptr != requiredAddr) {
        throw hal_exception("unable to remap file at same address");
    }
    return ptr;
}

/* unmap file, if mapped */
void hal::MMapFileLocal::unmapFile() {
    if (_basePtr != NULL) {
        if (::munmap(const_cast<void *>(_basePtr), _fileSize) < 0) {
            throw hal_errno_exception(_alignmentPath, "munmap failed", errno);
        }
        _basePtr = NULL;
    }
}

/* open the file for read access */
void hal::MMapFileLocal::openRead() {
    _fd = openFile();
    _fileSize = getFileStatSize(_fd);
    _basePtr = mapFile();
    loadHeader(false);
}

/* open the file for write access */
void hal::MMapFileLocal::openWrite(size_t fileSize) {
    _fd = openFile();
    if (_mode & CREATE_ACCESS) {
        adjustFileSize(0); // clear out existing data
        adjustFileSize(fileSize);
    } else if (_mode & WRITE_ACCESS) {
        adjustFileSize(getFileStatSize(_fd) + fileSize);
    }
    _basePtr = mapFile();
    if (_mode & CREATE_ACCESS) {
        createHeader();
    } else {
        loadHeader(true);
    }
}

/* close the file if open */
void hal::MMapFileLocal::closeFile() {
    if (_fd >= 0) {
    try_close:
        if (::close(_fd) < 0) {
            throw hal_errno_exception(_alignmentPath, "close failed", errno);
        }
        _fd = -1;
    }
}

#ifdef ENABLE_UDC
namespace hal {
    /* Class that implements UDC file version of MMapFile */
    class MMapFileUdc : public MMapFile {
      public:
        MMapFileUdc(const std::string &alignmentPath, unsigned mode, size_t fileSize);
        virtual void close();
        virtual ~MMapFileUdc();
        virtual bool isUdcProtocol() const {
            return true;
        }

      protected:
        virtual void fetch(size_t offset, size_t accessSize) const;

      private:
        struct udc2File *_udcFile;
    };
}

/* Constructor. Open or create the specified file. */
hal::MMapFileUdc::MMapFileUdc(const std::string &alignmentPath, unsigned mode, size_t fileSize)
    : MMapFile(alignmentPath, mode, true), _udcFile(NULL) {
    if (_mode & WRITE_ACCESS) {
        throw hal_exception("write access not supported for UDC:" + alignmentPath);
    }
    _udcFile = udc2FileMayOpen(const_cast<char *>(alignmentPath.c_str()), NULL, UDC_BLOCK_SIZE);
    if (_udcFile == NULL) {
        throw hal_exception("can't open " + alignmentPath);
    }
    udc2MMap(_udcFile);

    // get base point and fetch header
    _basePtr = udc2MMapFetch(_udcFile, 0, sizeof(MMapHeader));
    _fileSize = udc2SizeFromCache(const_cast<char *>(_alignmentPath.c_str()), NULL);
    loadHeader(false);
}

/* close file, marking as clean.  Don't  */
void hal::MMapFileUdc::close() {
    if (_basePtr == NULL) {
        throw hal_exception(_alignmentPath + ": MMapFile::close() called on closed file");
    }
    udc2FileClose(&_udcFile);
}

/* Destructor. write fields to header and close.  If write access and close
 * has not been called, file will me left mark dirty */
hal::MMapFileUdc::~MMapFileUdc() {
    if (_udcFile != NULL) {
        udc2FileClose(&_udcFile);
    }
}

/* fetch into UDC cache */
void hal::MMapFileUdc::fetch(size_t offset, size_t accessSize) const {
    if ((offset < _fileSize) and (offset + accessSize) > _fileSize) {
        // FIXME  - length off end, iterator does this
        accessSize = _fileSize - offset;
    }

    udc2MMapFetch(_udcFile, offset, accessSize);
}

#endif

/** create a MMapFile object, opening a local file */
hal::MMapFile *hal::MMapFile::factory(const std::string &alignmentPath, unsigned mode, size_t fileSize) {
    if (isUrl(alignmentPath)) {
        if (mode & (CREATE_ACCESS | WRITE_ACCESS)) {
            throw hal_exception("create or write access not support with URL: " + alignmentPath);
        }
#ifdef ENABLE_UDC
        return new MMapFileUdc(alignmentPath, mode, fileSize);
#else
        throw hal_exception("URL access requires UDC support to be compiled into HAL library: " + alignmentPath);
#endif
    } else {
        return new MMapFileLocal(alignmentPath, mode, fileSize);
    }
}
