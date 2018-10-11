#ifndef _MMAPFILE_H
#include "halDefs.h"
#include <cstddef>
#include <string>
#include <assert.h>
#include "halAlignmentInstance.h"

namespace hal {
    /* header for the file */
    struct mmapHeader {
        char format[32];
        char version[32];
        size_t nextOffset;
        size_t rootOffset;
        bool dirty;
    };
    typedef struct mmapHeader mmapHeader;

    /**
     * An mmapped HAL file.  This handles creation and opening of mapped
     * file.  
     * WARNING: When writing, close() must be explicitly called or file will
     * be left marked as dirty.
     */
    class MmapFile {
        public:

       
        const std::string& getStorageFormat() const {
            return STORAGE_FORMAT_MMAP;
        }

        inline size_t getRootOffset() const;
        inline void *toPtr(size_t offset,
                           size_t accessSize);
        inline const void *toPtr(size_t offset,
                                 size_t accessSize) const;
        inline size_t allocMem(size_t size,
                        bool isRoot=false);
        virtual ~MmapFile() {
        }
        
        protected:
        MmapFile(const std::string alignmentPath,
                 unsigned mode);
        /** close marks as clean, don't call on error, just delete */
        virtual void close() = 0;
        virtual void fetch(size_t offset,
                           size_t accessSize) const = 0;

        inline size_t alignRound(size_t size) const;
        void setHeaderPtr();
        void createHeader();
        void loadHeader(bool markDirty);
        void validateWriteAccess() const;
        void growFile(size_t size);
        virtual void growFileImpl(size_t size);
        inline void fetchIfNeeded(size_t offset,
                                  size_t accessSize) const;
        
        const std::string _alignmentPath;   // name of file for errors
        unsigned _mode;       // access mode
        void *_basePtr;       // location file is mapped
        mmapHeader *_header;  // pointer to header
        size_t _fileSize;     // size of file
        bool _mustFetch;      // fetch must be called on each access.

        private:
        MmapFile() {
            // no copying
        }

        static MmapFile* localFactory(const std::string& alignmentPath,
                                      unsigned mode = HAL_READ,
                                      size_t initSize = MMAP_DEFAULT_INIT_SIZE,
                                      size_t growSize = MMAP_DEFAULT_GROW_SIZE);
    };
}

/** Get the offset of the root object */
size_t hal::MmapFile::getRootOffset() const {
    assert(_header->rootOffset > 0);
    return _header->rootOffset;
}

/* fetch the range if required, else inline no-op */
void hal::MmapFile::fetchIfNeeded(size_t offset,
                                  size_t accessSize) const {
    if (_mustFetch) {
        fetch(offset, accessSize);
    }
}

/** Get pointer to the root a pointer.  Where accessSize is the
 * number of bytes that will be accessed, which is used when
 * pre-fetching is needed. If accessing an array, accessSize is size
 * of element, not the entire array.*/
void *hal::MmapFile::toPtr(size_t offset,
                           size_t accessSize) {
    fetchIfNeeded(offset, accessSize);
    return static_cast<char*>(_basePtr) + offset;
}

const void *hal::MmapFile::toPtr(size_t offset,
                                 size_t accessSize) const {
    fetchIfNeeded(offset, accessSize);
    return static_cast<const char*>(_basePtr) + offset;
}

/** Allocate new memory, resize file if necessary. If isRoot
 * is specified, it is stored as the root.  If AUTO_GROW mode,
 * then any pointers previously returned may become invalid. Thus
 * only offsets can be stored. */
size_t hal::MmapFile::allocMem(size_t size,
                               bool isRoot) {
    validateWriteAccess();
    if (_header->nextOffset + size > _fileSize) {
        growFile(size);
    }
    size_t offset = _header->nextOffset;
    _header->nextOffset += alignRound(size);
    if (isRoot) {
        _header->rootOffset = offset;
    }
    return offset;
}

/* round up to alignment size */
size_t hal::MmapFile::alignRound(size_t size) const {
    return ((size + (sizeof(size_t) - 1)) / sizeof(size_t)) * sizeof(size_t);
}
#endif

// Local Variables:
// mode: C ++
// End
