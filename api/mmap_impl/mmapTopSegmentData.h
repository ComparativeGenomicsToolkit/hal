#ifndef _MMAPTOPSEGMENTDATA_H
#define _MMAPTOPSEGMENTDATA_H

namespace hal {
class MMapTopSegmentData
{
    public:
    void setStartPosition(hal_index_t startPosition) { _startPosition = startPosition; };
    void setBottomParseIndex(hal_index_t parseIndex) { _bottomParseIndex = parseIndex; };
    void setNextParalogyIndex(hal_index_t paralogyIndex) { _paralogyIndex = paralogyIndex; };
    void setParentIndex(hal_index_t parentIndex) { _parentIndex = parentIndex; };
    void setReversed(bool reversed) { _reversed = reversed; };

    hal_index_t getStartPosition() const { return _startPosition; };
    hal_index_t getBottomParseIndex() const { return _bottomParseIndex; };
    hal_index_t getNextParalogyIndex() const { return _paralogyIndex; };
    hal_index_t getParentIndex() const { return _parentIndex; };
    hal_index_t getReversed() const { return _reversed; };

    private:
    hal_index_t _startPosition;
    hal_index_t _bottomParseIndex;
    hal_index_t _paralogyIndex;
    hal_index_t _parentIndex;
    bool _reversed;
};
}
#endif
// Local Variables:
// mode: c++
// End:
