#ifndef _MMAPSEQUENCEDATA_H
#define _MMAPSEQUENCEDATA_H
#include "mmapAlignment.h"
namespace hal {
class MMapSequenceData {
    friend class MMapSequence;
    public:
    const char *getName(MMapAlignment *alignment) const { return (const char *) alignment->resolveOffset(_nameOffset, _nameLength); };
    void setName(MMapAlignment *alignment, const std::string &newName) {
        size_t size = newName.size() + 1;
        _nameOffset = alignment->allocateNewArray(sizeof(char) * size);
        strncpy((char *) alignment->resolveOffset(_nameOffset, size), newName.c_str(), size);
        _nameLength = size;
    };
    protected:
    hal_index_t _startPosition;
    hal_index_t _index;
    hal_size_t _length;
    hal_index_t _topSegmentStartIndex;
    hal_index_t _bottomSegmentStartIndex;
    hal_size_t _numTopSegments;
    hal_size_t _numBottomSegments;
    size_t _nameLength;
    size_t _nameOffset;
};
}
#endif
