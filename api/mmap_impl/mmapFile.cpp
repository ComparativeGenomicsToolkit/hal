#include "mmapFile.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>


/* constants for header */
static const std::string FORMAT_NAME = "MMAP";
static const std::string MMAP_VERSION = "1.0";


/* get the file size from the OS */
static size_t getFileStatSize(int fd) {
    struct stat fileStat;
    if (::fstat(fd, &fileStat) < 0) {
        throw hal_errno_exception("stat failed", errno);
    }
    return fileStat.st_size;
}

/* constructor, used only by derived classes */
hal::MmapFile::MmapFile(const std::string alignmentPath,
                        unsigned mode):
    _alignmentPath(alignmentPath),  _mode(halDefaultAccessMode(mode)),
    _basePtr(NULL), _fileSize(0), _mustFetch(false) {
}

/* error if file is not open for write accecss */
void hal::MmapFile::validateWriteAccess() const {
    if ((_mode & WRITE_ACCESS) == 0) {
        throw hal_exception(_alignmentPath + " is not open for write access");
    }
}

/* setup pointer to header */
void hal::MmapFile::setHeaderPtr() {
    fetchIfNeeded(0, sizeof(mmapHeader));
    _header = static_cast<mmapHeader*>(_basePtr);
}

/* validate the file header and save a pointer to it. */
void hal::MmapFile::loadHeader(bool markDirty) {
    if (_fileSize < sizeof(mmapHeader)) {
        throw hal_exception(_alignmentPath + ": file size of " + std::to_string(_fileSize)
                            + " is less that header size of " + std::to_string(sizeof(mmapHeader)));
    }
    setHeaderPtr();

    // don't print found strings as it might have garbage
    if (::strcmp(_header->format, FORMAT_NAME.c_str()) != 0) {
        throw hal_exception(_alignmentPath + ": invalid file header, expected format name of '" +
                            FORMAT_NAME + "'");
    }
    // FIXME: to need to check version compatibility
    if (::strcmp(_header->mmapVersion, MMAP_VERSION.c_str()) != 0) {
        throw hal_exception(_alignmentPath + ": incompatible mmap format versions: "
                            + "file version " + _header->mmapVersion
                            + ", mmap API version " + MMAP_VERSION);
    }
    if ((_header->nextOffset < sizeof(mmapHeader)) or (_header->nextOffset > _fileSize)) {
        throw hal_exception(_alignmentPath + ": header nextOffset field out of bounds, probably file corruption");
    }
    if ((_header->rootOffset < sizeof(mmapHeader)) or (_header->rootOffset > _fileSize)) {
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
void hal::MmapFile::createHeader() {
    assert(_mode & WRITE_ACCESS);
    setHeaderPtr();
    assert(strlen(FORMAT_NAME) < sizeof(_header->format));
    strncpy(_header->format, FORMAT_NAME.c_str(), sizeof(_header->format)-1);
    assert(strlen(MMAP_VERSION) < sizeof(_header->version));
    strncpy(_header->mmapVersion, MMAP_VERSION.c_str(), sizeof(_header->mmapVersion)-1);
    assert(strlen(HAL_VERSION) < sizeof(_header->halVersion));
    strncpy(_header->halVersion, HAL_VERSION.c_str(), sizeof(_header->halVersion)-1);
    _header->nextOffset = alignRound(sizeof(mmapHeader));
    _header->dirty = true;
    _header->nextOffset = _header->nextOffset;
}

/* Grow file to allow for at least the specified amount.  This remaps the *
 * file, so the same sufficient virtual space must be available at the address
 * and it is expensive. */
void hal::MmapFile::growFile(size_t size) {
    assert(_mode & WRITE_ACCESS);
    growFileImpl(size);
}

/* override this for classes that support growing file */
void hal::MmapFile::growFileImpl(size_t size) {
    throw hal_exception("logic error: growFile() not available for this MmapFile implementation");
}

namespace hal {
    /* Class that implements local file version of MmapFile */
    class MmapFileLocal: public MmapFile {
        public:
        MmapFileLocal(const std::string alignmentPath,
                      unsigned mode,
                      size_t initSize,
                      size_t growSize);
        virtual void close();
        virtual ~MmapFileLocal();

        protected:
        virtual void fetch(size_t offset,
                           size_t accessSize) const {
            // no-op
        }

        private:
        int openFile();
        void closeFile();
        void adjustFileSize(size_t size);
        void* mapFile(void *requiredAddr=NULL);
        void unmapFile();
        void openRead();
        void openWrite(size_t initSize);
        void growFileImpl(size_t size);

        int _fd;              // open file descriptor
        size_t _growSize;     // amount to grow file by when needed.
    };
}


/* Constructor. Open or create the specified file. */
hal::MmapFileLocal::MmapFileLocal(const std::string alignmentPath,
                                  unsigned mode,
                                  size_t initSize,
                                  size_t growSize):
    MmapFile(alignmentPath, mode), _fd(-1), _growSize(growSize) {
    if (_mode & WRITE_ACCESS) {
        openWrite(initSize);
    } else {
        openRead();
    }
}

/* close file, marking as clean.  Don't  */
void hal::MmapFileLocal::close() {
    if (_basePtr == NULL) {
        throw hal_exception(_alignmentPath + ": MmapFile::close() called on closed file");
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
hal::MmapFileLocal::~MmapFileLocal() {
    unmapFile();
    closeFile();
}

/* open the file for the specified mode */
int hal::MmapFileLocal::openFile() {
    assert(_fd < 0);
    unsigned openMode = 0;
    if (_mode & WRITE_ACCESS) {
        openMode = O_RDWR | ((_mode & CREATE_ACCESS) ? (O_CREAT|O_TRUNC) : 0);
    } else {
        openMode = O_RDONLY;
    }
    int fd = ::open(_alignmentPath.c_str(), openMode, 0777);
    if (fd < 0) {
        throw hal_errno_exception(_alignmentPath, "open failed", errno);
    }
    return fd;
}

/* change size size of the file, possibly deleting data. */
void hal::MmapFileLocal::adjustFileSize(size_t size) {
    if (ftruncate(_fd, size) < 0) {
        throw hal_errno_exception(_alignmentPath, "set size failed", errno);
    }
    _fileSize = size;
}

/* map file into memory */
void* hal::MmapFileLocal::mapFile(void *requiredAddr) {
    assert(_basePtr == NULL);
    unsigned prot = PROT_READ | ((_mode & WRITE_ACCESS) ? PROT_WRITE : 0);
    void *ptr = mmap(requiredAddr, _fileSize, prot, MAP_SHARED|MAP_FILE, _fd, 0);
    if (ptr == MAP_FAILED) {
        throw hal_errno_exception(_alignmentPath, "mmap failed", errno);
    }
    return ptr;
}

/* unmap file, if mapped */
void hal::MmapFileLocal::unmapFile() {
    if (_basePtr != NULL) {
        if (::munmap(const_cast<void*>(_basePtr), _fileSize) < 0) {
            throw hal_errno_exception(_alignmentPath, "munmap failed", errno);
        }
        _basePtr = NULL;
    }
}

/* open the file for read access */
void hal::MmapFileLocal::openRead() {
    _fd = openFile();
    _fileSize = getFileStatSize(_fd);
    _basePtr = mapFile();
    loadHeader(false);
}

/* open the file for write access */
void hal::MmapFileLocal::openWrite(size_t initSize) {
    _fd = openFile();
    if (_mode & CREATE_ACCESS) {
        adjustFileSize(0);  // clear out existing data
    }
    if (initSize > _fileSize) {
        adjustFileSize(initSize);
    }
    _basePtr = mapFile();
    if (_mode & CREATE_ACCESS) {
        createHeader();
    } else {
        loadHeader(true);
    }
}

/* close the file if open */
void hal::MmapFileLocal::closeFile() {
    if (_fd >= 0) {
        if (::close(_fd) < 0) {
            throw hal_errno_exception(_alignmentPath, "close failed", errno);
        }
        _fd = -1;
    }
}

/* grow file  */
void hal::MmapFileLocal::growFileImpl(size_t size) {
    assert(_mode & WRITE_ACCESS);
    size_t newSize = _growSize;
    if (newSize < size) {
        newSize += size;  // will leave in extra
    }
    void * requiredAddr = _basePtr;
    unmapFile();
    mapFile(requiredAddr);
}

/** create a MmapFile object, opening a local file */
hal::MmapFile *hal::MmapFile::localFactory(const std::string& alignmentPath,
                                           unsigned mode,
                                           size_t initSize,
                                           size_t growSize) {
    return new MmapFileLocal(alignmentPath, mode, initSize, growSize);
}

