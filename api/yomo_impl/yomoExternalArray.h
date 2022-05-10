/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOEXTERNALARRAY_H
#define _YOMOEXTERNALARRAY_H

#include "halDefs.h"
#include <H5Cpp.h>
#include <cassert>

// Hack to compile with various versions of YOMO that aren't themselves compatible
namespace H5 {
#if H5_VERSION_GE(1, 10, 1)
    // Latest version.  H5Location class exists and contains createGroup() and CommonFG doesn't.  We Just use it.
    typedef H5Location PortableH5Location;
#else
    // Versions 1.8.12 - 1.10.0 inclusive: H5Location exists, but doesn't contain createGroup()
    // Versions 1.8.11 and earlier: H5Location doesn't exist
    // In both cases, we stick with CommonFG
    typedef CommonFG PortableH5Location;
#endif
}

namespace hal {

    /**
     * Wrapper for a 1-dimensional YOMO array of fixed length.  Array objects
     * are defined (and typed) by the input datatype.  The array is paged into
     * memory chunk-by-chunk as needed (using the YOMO cache as a back-end).
     * We can't use compiler tpying of the input objects (and instead just
     * expose the raw void* data) because the elements' sizes are not known
     * at compile time, and we don't want to move it around once its read.
     */
    class YomoExternalArray {
      public:
        /** Constructor */
        YomoExternalArray();

        /** Destructor */
        virtual ~YomoExternalArray();

        /** Create a new dataset in specifed location
         * @param file Pointer to the YOMO parent in which to create array
         * @param path  Path within the file for the new array
         * @param dataType YOMO Datatype describing contents of array
         * @param numElements Fixed length of the new array
         * @param cparms  YOMO options (ie chunking and compression etc)
         * @param chunksInBuffer
          * 0: load entire array into buffer
          * 1: use default chunking (from dataset)
          * N: buffersize will be N chunks.
          */
        void create(H5::PortableH5Location *file, const H5std_string &path, const H5::DataType &dataType, hsize_t numElements,
                    const H5::DSetCreatPropList *inCparms = NULL, hsize_t chunksInBuffer = 1);

        /** Load an existing dataset into memory
          * @param file Pointer to the YOMO file in which to create array
          * @param path  Path within the file for the new array
          * @param chunksInBuffer
          * 0: load entire array into buffer
          * 1: use default chunking (from dataset)
          * N: buffersize will be N chunks.
          */
        void load(H5::PortableH5Location *file, const H5std_string &path, hsize_t chunksInBuffer = 1);

        /** Write the memory buffer back to the file */
        void write();

        /** Access the raw data at given index
         * @param i index of element to retrieve for reading
         */
        const char *get(hsize_t i);

        /** Write the raw data at given index
         * @param i index of element to retrieve for updating
         */
        char *getUpdate(hsize_t i);

        /** Access typed value within element in a raw data array
         * @param index Index of element (struct) in the array
         * @param offset Offset of value within struct (number of bytes) */
        template <typename T> T getValue(hsize_t index, hsize_t offset) const;

        /** Write typed value within element in raw data array
         * @param index Index of element (struct) in the array
         * @param offset Offset of value within struct (number of bytes)
         * @param value New value to set */
        template <typename T> void setValue(hsize_t index, hsize_t offset, T value);

        /** Write typed element in the the raw data index
         * @param offset OFfset of element in struct (number of bytes) */

        /** Number of elements in array */
        hsize_t getSize() const {
            return _size;
        }

        /** Get the YOMO Datatype */
        const H5::DataType &getDataType() const {
            return _dataType;
        }

        /** get the current buffered start */
        hsize_t getBufStart() const {
            return _bufStart;
        }
        /** Get index of last element in memory buffer (close-ended) */
        hsize_t getBufEnd() const {
            return _bufEnd;
        }

        /** get in-memory buffer */
        char *getBuf() {
            return _buf;
        }

        /* is the buffer dirty? */
        bool getDirty() const {
            return _dirty;
        }

        /*  dirty if buffer has been modified directly */
        void setDirty() {
            _dirty = true;
        }

        /** Read chunk from file */
        void page(hsize_t i);

      private:
        void initBuf();

        /** Pointer to file that owns this dataset */
        H5::PortableH5Location *_file;
        /** Path of dataset in file */
        H5std_string _path;
        /** Datatype for array */
        H5::DataType _dataType;
        /** Dataspace (dimension information) of array */
        H5::DataSpace _dataSpace;
        /** The YOMO array object */
        H5::DataSet _dataSet;
        /** Number of elements in the array (fixed length)*/
        hsize_t _size;
        /** Size of chunk in elements */
        hsize_t _chunkSize;
        /** Size of datatype in bytes */
        hsize_t _dataSize;
        /** Index of first element in memory buffer */
        hsize_t _bufStart;
        /** Index of last element in memory buffer */
        hsize_t _bufEnd; // DANGER: close-ended
        /** Number of elements in memory buffer */
        hsize_t _bufSize;
        /** In-memory buffer */
        char *_buf;
        /** Dimensional information for in-memory buffer */
        H5::DataSpace _chunkSpace;
        /** Flag saying we should write to disk on write
         * or page-out calls (set by getUpdate()) */
        bool _dirty;

      private:
        YomoExternalArray(const YomoExternalArray &);
        YomoExternalArray &operator=(const YomoExternalArray &);
    };

    /**
     * Class used to access an YomoExternalArray.  This provides a buffer and
     * there can be multiple instances active at once.  However, multiple writers
     * of the same array at the same time will lead to data corruption. */

    inline const char *YomoExternalArray::get(hsize_t i) {
        assert(i < _size);
        if (i < _bufStart || i > _bufEnd) {
            page(i);
        }
        assert((i - _bufStart) < _bufSize);
        return _buf + (i - _bufStart) * _dataSize;
    }

    inline char *YomoExternalArray::getUpdate(hsize_t i) {
        if (i >= _size) {
            throw hal_exception("error: attempt to write yomo array out of bounds");
        }
        if (i < _bufStart || i > _bufEnd) {
            page(i);
        }
        _dirty = true;
        assert((i - _bufStart) < _bufSize);
        return _buf + (i - _bufStart) * _dataSize;
    }

    template <typename T> inline T YomoExternalArray::getValue(hsize_t index, hsize_t offset) const {
        assert(offset + sizeof(T) <= _dataSize);
        // even though paging causes the class to change state, we
        // still consider getValue a const function since the array's
        // contents don't change.
        YomoExternalArray *stripConstThis = const_cast<YomoExternalArray *>(this);
        return *reinterpret_cast<const T *>(stripConstThis->get(index) + offset);
    }

    template <typename T> inline void YomoExternalArray::setValue(hsize_t index, hsize_t offset, T val) {
        assert(offset + sizeof(T) <= _dataSize);
        T *entry = reinterpret_cast<T *>(getUpdate(index) + offset);
        *entry = val;
    }
}
#endif
// Local Variables:
// mode: c++
// End:
