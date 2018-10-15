#ifndef _MMAPBOTTOMSEGMENTDATA_H
#define _MMAPBOTTOMSEGMENTDATA_H

namespace hal {
class MMapBottomSegmentData
{
    public:
    void setStartPosition(hal_index_t startPosition) { _startPosition = startPosition; };
    void setTopParseIndex(hal_index_t parseIndex) { _topParseIndex = parseIndex; };
    void setChildIndex(hal_size_t child, hal_index_t childIndex) {
        *getChildIndexLocation(child) = childIndex;
    };
    void setChildReversed(hal_size_t numChildren, hal_size_t child, bool childReversed) {
        *getChildReversedLocation(numChildren, child) = childReversed;
    };

    hal_index_t getStartPosition() const { return _startPosition; };
    hal_index_t getTopParseIndex() const { return _topParseIndex; };
    hal_index_t getChildIndex(hal_size_t child) const {
        return *getChildIndexLocation(child);
    };
    hal_index_t getChildReversed(hal_size_t numChildren, hal_size_t child) const {
        return *getChildReversedLocation(numChildren, child);
    };

    // Get on-disk size of this element for the given genome. NB: the size is
    // rounded up to the next 8-byte boundary for alignment purposes.
    static size_t getSize(const Genome *genome) {
        size_t extraAlignmentBytes = 0;
        if ((genome->getNumChildren() % 8) != 0) {
            extraAlignmentBytes = 8 - (genome->getNumChildren() % 8);
        }
        return sizeof(hal_index_t) * (2 + genome->getNumChildren()) + (genome->getNumChildren() + extraAlignmentBytes);
    };

    private:
    hal_index_t *getChildIndexLocation(hal_size_t child) const {
        return const_cast<hal_index_t *>(&_topParseIndex + child);
    }
    bool *getChildReversedLocation(hal_size_t numChildren, hal_size_t child) const {
        return ((bool *) getChildIndexLocation(numChildren)) + child;
    }
    hal_index_t _startPosition;
    hal_index_t _topParseIndex;
};
}
#endif
