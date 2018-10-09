#include "mmapFile.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/mman.h>


/* constants for header */
static const char *FORMAT_NAME = "MMAP";
static const char *FORMAT_VERSION = "1.0";


/* get the file size from the OS */
static size_t getFileStatSize(int fd) {
    struct stat fileStat;
    if (::fstat(fd, &fileStat) < 0) {
        throw hal_errno_exception("stat failed", errno);
    }
    return fileStat.st_size;
}

/* constructor, used only by derived classes */
hal::MmapFile::MmapFile(const std::string fileName,
                        unsigned mode):
    _fileName(fileName),  _mode(mode), _basePtr(NULL), _fileSize(0), _mustFetch(false) {
    // make mode sane and validate
    if (_mode & MMAP_CREATE) {
        _mode |= MMAP_WRITE;
    }
    if (_mode & MMAP_WRITE) {
        _mode |= MMAP_READ;
    }
    if ((_mode & (MMAP_READ|MMAP_WRITE|MMAP_CREATE)) == 0) {
        throw hal_exception(fileName + ": must specify at least one of READ, WRITE, or CREATE on open");
    }
}

/* error if file is not open for write accecss */
void hal::MmapFile::validateWriteAccess() const {
    if ((_mode & MMAP_WRITE) == 0) {
        throw hal_exception(_fileName + " is not open for write access");
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
        throw hal_exception(_fileName + ": file size of " + std::to_string(_fileSize)
                            + " is less that header size of " + std::to_string(sizeof(mmapHeader)));
    }
    setHeaderPtr();

    // don;t print found strings as it might have garbage
    if (::strcmp(_header->format, FORMAT_NAME) != 0) {
        throw hal_exception(_fileName + ": invalid file header, expected format name of '" +
                            FORMAT_NAME + "'");
    }
    if (::strcmp(_header->version, FORMAT_VERSION) != 0) {
        throw hal_exception(_fileName + ": invalid file header, expected version of '" +
                            FORMAT_VERSION + "'");
    }
    if ((_header->nextOffset < sizeof(mmapHeader)) or (_header->nextOffset > _fileSize)) {
        throw hal_exception(_fileName + ": header nextOffset field out of bounds, probably file corruption");
    }
    if ((_header->rootOffset < sizeof(mmapHeader)) or (_header->rootOffset > _fileSize)) {
        throw hal_exception(_fileName + ": header rootOffset field out of bounds, probably file corruption");
    }
    if (_header->dirty) {
        throw hal_exception(_fileName + ": file is marked as dirty, most likely an inconsistent state.");
    }
    if (markDirty) {
        _header->dirty = true;
    }
}

/* create the header */
void hal::MmapFile::createHeader() {
    assert(_mode & MMAP_WRITE);
    setHeaderPtr();
    assert(strlen(FORMAT_NAME) < sizeof(_header->format));
    strncpy(_header->format, FORMAT_NAME, sizeof(_header->format));
    assert(strlen(FORMAT_VERSION) < sizeof(_header->version));
    strncpy(_header->version, FORMAT_VERSION, sizeof(_header->version));
    _header->nextOffset = alignRound(sizeof(mmapHeader));
    _header->dirty = true;
    _header->nextOffset = _header->nextOffset;
}

/* Grow file to allow for at least the specified amount.  This remaps the
 * file, so it is expensive and existing pointer become invalid. */
void hal::MmapFile::growFile(size_t size) {
    assert(_mode & MMAP_WRITE);
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
        MmapFileLocal(const std::string fileName,
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
        void* mapFile();
        void unmapFile();
        void openRead();
        void openWrite(size_t initSize);
        void growFileImpl(size_t size);

        int _fd;              // open file descriptor
        size_t _growSize;     // amount to grow file by when needed.
    };
}


/* Constructor. Open or create the specified file. */
hal::MmapFileLocal::MmapFileLocal(const std::string fileName,
                                  unsigned mode,
                                  size_t initSize,
                                  size_t growSize):
    MmapFile(fileName, mode), _fd(-1), _growSize(growSize) {
    if (_mode & MMAP_WRITE) {
        openWrite(initSize);
    } else {
        openRead();
    }
}

/* close file, marking as clean.  Don't  */
void hal::MmapFileLocal::close() {
    if (_basePtr == NULL) {
        throw hal_exception(_fileName + ": MmapFile::close() called on closed file");
    }
    if (_mode & MMAP_WRITE) {
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
    if (_mode & MMAP_WRITE) {
        openMode = O_RDWR | ((_mode & MMAP_CREATE) ? (O_CREAT|O_TRUNC) : 0);
    } else {
        openMode = O_RDONLY;
    }
    int fd = ::open(_fileName.c_str(), openMode);
    if (fd < 0) {
        throw hal_errno_exception(_fileName, "open failed", errno);
    }
    return fd;
}

/* change size size of the file, possibly deleting data. */
void hal::MmapFileLocal::adjustFileSize(size_t size) {
    if (ftruncate(_fd, size) < 0) {
        throw hal_errno_exception(_fileName, "set size failed", errno);
    }
    _fileSize = size;
}

/* map file into memory */
void* hal::MmapFileLocal::mapFile() {
    assert(_basePtr == NULL);
    unsigned prot = PROT_READ | ((_mode & MMAP_WRITE) ? PROT_WRITE : 0);
    void *ptr = mmap(0, _fileSize, prot, MAP_SHARED|MAP_FILE, _fd, 0);
    if (ptr == MAP_FAILED) {
        throw hal_errno_exception(_fileName, "mmap failed", errno);
    }
    return ptr;
}

/* unmap file, if mapped */
void hal::MmapFileLocal::unmapFile() {
    if (_basePtr != NULL) {
        if (::munmap(const_cast<void*>(_basePtr), _fileSize) < 0) {
            throw hal_errno_exception(_fileName, "munmap failed", errno);
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
    if (_mode & MMAP_CREATE) {
        adjustFileSize(0);  // clear out existing data
    }
    if (initSize > _fileSize) {
        adjustFileSize(initSize);
    }
    _basePtr = mapFile();
    if (_mode & MMAP_CREATE) {
        createHeader();
    } else {
        loadHeader(true);
    }
}

/* close the file if open */
void hal::MmapFileLocal::closeFile() {
    if (_fd >= 0) {
        if (::close(_fd) < 0) {
            throw hal_errno_exception(_fileName, "close failed", errno);
        }
        _fd = -1;
    }
}

/* grow file  */
void hal::MmapFileLocal::growFileImpl(size_t size) {
    size_t newSize = _growSize;
    if (newSize < size) {
        newSize += size;  // will leave in extra
    }
    
}

/** create a MmapFile object, opening a local file */
hal::MmapFile* hal::MmapFile::localFactory(const std::string& fileName,
                                           unsigned mode,
                                           size_t initSize,
                                           size_t growSize) {
    return new MmapFileLocal(fileName, mode, initSize, growSize);
}

