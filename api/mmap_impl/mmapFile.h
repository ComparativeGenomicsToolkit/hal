#ifndef _MMAPFILE_H
#include <cstddef>
#include <string>

namespace hal
{
struct mmapHeader;
typedef struct mmapHeader mmapHeader;

/**
 * An mmapped HAL file.  This handles creation and opening of mapped
 * file.
 */
class MmapFile
{
    public:
    // Open modes
    enum {
        READ = 0x01,      // read-access
        WRITE = 0x02,     // write-access
        CREATE = 0x04,    // initialize a new file, must be empty if it exists
        AUTO_GROW = 0x08  // allow auto-growing, see warnings
    };

    off_t getRootOffset() const;
    void *toPtr(off_t offset,
                size_t accessSize);
    const void *toPtr(off_t offset,
                      size_t accessSize) const;
    size_t allocMem(size_t size,
                    bool isRoot=false);

private:
    size_t alignRound(size_t size) const;
    mmapHeader* getHeader();
    const mmapHeader* getHeader() const;
    void validateWriteAccess() const;
    void growFile(size_t size);
    void storeRoot(off_t rootOffset);

    const std::string _fileName;   // name of file for errors
    bool _writeAccess;    // is it open for write access */
    int _fd;              // open file descriptor
    void *_basePtr;       // location file is mapped
    size_t _fileSize;     // size of file
    size_t _growSize;     // size to grow file on resize
    off_t _rootOffset;    // offset of root object
    off_t _nextOffset;    // next byte to allocated when creating (is aligned)
};


/** Get the offset of the root object */
off_t MmapFile::getRootOffset() const
{
    return _rootOffset;
}

/** Get pointer to the root a pointer.  Where accessSize is the
 * number of bytes that will be accessed, which is used when
 * pre-fetching is needed. If accessing an array, accessSize is size
 * of element, not the entire array.*/
void *MmapFile::toPtr(off_t offset,
                      size_t accessSize)
{
    return static_cast<char*>(_basePtr) + offset;
}

const void *MmapFile::toPtr(off_t offset,
                            size_t accessSize) const
{
    return static_cast<const char*>(_basePtr) + offset;
}

/** Allocate new memory, resize file if necessary. If isRoot
 * is specified, it is stored as the root.  If AUTO_GROW mode,
 * then any pointers previously returned may become invalid. Thus
 * only offsets can be stored. */
size_t MmapFile::allocMem(size_t size,
                          bool isRoot)
{
    validateWriteAccess();
    if (_nextOffset + size > _fileSize)
    {
        growFile(size);
    }
    size_t offset = _nextOffset;
    _nextOffset += alignRound(size);
    if (isRoot)
    {
        storeRoot(offset);
    }
    return offset;
}

/* round up to alignment size */
size_t MmapFile::alignRound(size_t size) const
{
    return ((size + (sizeof(size_t) - 1)) / sizeof(size_t)) * sizeof(size_t);
}
}
#endif

// Local Variables:
// mode: C ++
// End
