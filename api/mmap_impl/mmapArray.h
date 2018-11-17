#ifndef _MMAPARRAY_H
#define _MMAPARRAY_H
#include "halMetaData.h"

namespace hal {
class MMapArrayData {
public:
    size_t _elementSize;
    size_t _capacity;
    size_t _length;
};

template <class T>
class MMapArray {
    public:
    // Construct & initalize new MMapArray.
    MMapArray(MMapAlignment *alignment) : _alignment(alignment), _data(NULL) { grow(8); };
    // Construct a MMapArray representing the existing array at this offset.
    MMapArray(MMapAlignment *alignment, size_t offset) : _alignment(alignment), _offset(offset) { _data = (MMapArrayData *) _alignment->resolveOffset(_offset, sizeof(MMapArrayData)); };
    T *operator[](size_t index) { return getSlice(index, 1); };
    T *getSlice(size_t index, size_t length) { return (T *) _alignment->resolveOffset(_offset + sizeof(MMapArrayData) + index * _data->_elementSize, _data->_elementSize * length); };
    size_t getOffset() { return _offset; };
    size_t getCapacity() { return _data->_capacity; };
    size_t getLength() { return _data->_length; };

    void grow(size_t capacity) {
        size_t size = sizeof(MMapArrayData) + capacity * sizeof(T);
        size_t newOffset = _alignment->allocateNewArray(size);
        MMapArrayData *newData = (MMapArrayData *) _alignment->resolveOffset(newOffset, size);
        if (_data != NULL) {
            memcpy(newData, _data, sizeof(MMapArrayData) + _data->_length * _data->_elementSize);
        }
        _data = newData;
        _offset = newOffset;
        _data->_capacity = capacity;
        _data->_elementSize = sizeof(T);
    };

    size_t setLength(size_t length) {
        if (length <= _data->_capacity) {
            _data->_length = length;
        } else {
            grow(length);
            _data->_length = length;
        }
        return getOffset();
    }
    private:
    MMapAlignment *_alignment;
    size_t _offset;
    MMapArrayData *_data;
};
}
#endif
// Local Variables:
// mode: c++
// End:
