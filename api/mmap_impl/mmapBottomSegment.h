#ifndef _MMAPBOTTOMSEGMENT_H
#define _MMAPBOTTOMSEGMENT_H
#include "halBottomSegment.h"
#include "halGenome.h"
#include "mmapGenome.h"
#include "mmapBottomSegmentData.h"
#include <cassert>

namespace hal {
class MMapBottomSegment : public BottomSegment
{
    public:
    MMapBottomSegment(const MMapGenome *genome, hal_index_t arrayIndex) :
        _index(arrayIndex) {
            _genome = const_cast<MMapGenome *>(genome);
            _data = _genome->getBottomSegmentPointer(arrayIndex);
        };

    // SEGMENT INTERFACE
    void setArrayIndex(Genome *genome, hal_index_t arrayIndex) {
        _genome = dynamic_cast<MMapGenome *>(genome);
        _data = _genome->getBottomSegmentPointer(arrayIndex);
        _index = arrayIndex;
    };
    void setArrayIndex(const Genome *genome, hal_index_t arrayIndex) const {
        const MMapGenome *mGenome = dynamic_cast<const MMapGenome *>(genome);
        _genome = const_cast<MMapGenome *>(mGenome);
        _data = _genome->getBottomSegmentPointer(arrayIndex);
        _index = arrayIndex;
    };
    const Genome* getGenome() const;
    Genome* getGenome();
    const Sequence* getSequence() const;
    hal_index_t getStartPosition() const { return _data->getStartPosition(); };
    hal_index_t getEndPosition() const;
    hal_size_t getLength() const;
    void getString(std::string& outString) const;
    void setCoordinates(hal_index_t startPos, hal_size_t length);
    hal_index_t getArrayIndex() const;
    bool leftOf(hal_index_t genomePos) const;
    bool rightOf(hal_index_t genomePos) const;
    bool overlaps(hal_index_t genomePos) const;
    bool isFirst() const;
    bool isLast() const;
    bool isMissingData(double nThreshold) const;
    bool isTop() const;
    hal_size_t getMappedSegments(
        MappedSegmentSet& outSegments,
        const Genome* tgtGenome,
        const std::set<const Genome*>* genomesOnPath,
        bool doDupes,
        hal_size_t minLength,
        const Genome *coalescenceLimit,
        const Genome *mrca) const;
    void print(std::ostream& os) const;

    // BOTTOM SEGMENT INTERFACE
    hal_size_t getNumChildren() const;
    hal_index_t getChildIndex(hal_size_t i) const { return _data->getChildIndex(i); };
    hal_index_t getChildIndexG(const Genome* childGenome) const;
    bool hasChild(hal_size_t child) const;
    bool hasChildG(const Genome* childGenome) const;
    void setChildIndex(hal_size_t i, hal_index_t childIndex) { _data->setChildIndex(i, childIndex); };
    bool getChildReversed(hal_size_t i) const { return _data->getChildReversed(_genome->getNumChildren(), i); };
    void setChildReversed(hal_size_t child, bool isReversed) { _data->setChildReversed(_genome->getNumChildren(), child, isReversed); };
    hal_index_t getTopParseIndex() const { return _data->getTopParseIndex(); };
    void setTopParseIndex(hal_index_t parseIndex) { _data->setTopParseIndex(parseIndex); };
    hal_offset_t getTopParseOffset() const;
    bool hasParseUp() const;
    hal_index_t getLeftChildIndex(hal_size_t i) const;
    hal_index_t getRightChildIndex(hal_size_t i) const;

    private:
    // Return a pointer to the data for the segment *after* this one in the array.
    MMapBottomSegmentData *getNextData() const {
        return (MMapBottomSegmentData *) (((char *) _data) + MMapBottomSegmentData::getSize(_genome));
    };
    mutable MMapBottomSegmentData *_data;
    mutable MMapGenome *_genome;
    mutable hal_index_t _index;
};

inline const Genome* MMapBottomSegment::getGenome() const
{
  return _genome;
}

inline Genome* MMapBottomSegment::getGenome()
{
  return _genome;
}

inline hal_index_t MMapBottomSegment::getEndPosition() const
{
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t MMapBottomSegment::getLength() const
{
  return getNextData()->getStartPosition() - _data->getStartPosition();
}

inline const Sequence* MMapBottomSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline hal_size_t MMapBottomSegment::getNumChildren() const
{
  return _genome->getNumChildren();
}

inline 
hal_index_t MMapBottomSegment::getChildIndexG(const Genome* childGenome) const
{
  assert(_index >= 0);
  return getChildIndex(_genome->getChildIndex(childGenome));
}

inline bool MMapBottomSegment::hasParseUp() const
{
  return getTopParseIndex() != NULL_INDEX;
}

inline bool MMapBottomSegment::hasChild(hal_size_t i) const
{
  return getChildIndex(i) != NULL_INDEX;
}

inline bool MMapBottomSegment::hasChildG(const Genome* childGenome) const
{
  return getChildIndexG(childGenome) != NULL_INDEX;
}

inline hal_index_t MMapBottomSegment::getArrayIndex() const
{
  return _index;
}

inline bool MMapBottomSegment::leftOf(hal_index_t genomePos) const
{
  return getEndPosition() < genomePos;
}

inline bool MMapBottomSegment::rightOf(hal_index_t genomePos) const
{
  return getStartPosition() > genomePos;
}

inline bool MMapBottomSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

inline bool MMapBottomSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getBottomSegmentArrayIndex();
}

inline bool MMapBottomSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == _genome->getNumBottomSegments() || 
     _index == getSequence()->getBottomSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumBottomSegments() - 1;
}

inline bool MMapBottomSegment::isTop() const
{
  return false;
}

inline hal_size_t MMapBottomSegment::getMappedSegments(
  MappedSegmentSet& outSegments,
  const Genome* tgtGenome,
  const std::set<const Genome*>* genomesOnPath,
  bool doDupes,
  hal_size_t minLength,
  const Genome *coalescenceLimit,
  const Genome *mrca) const
{
  throw hal_exception("Internal error. MMapSegment interface should "
                      "at some point go through the sliced segment");
}

inline hal_index_t MMapBottomSegment::getLeftChildIndex(hal_size_t i) const
{
  assert(isFirst() == false);
  MMapBottomSegment leftSeg(_genome, _index - 1);
  return leftSeg.getChildIndex(i);
}

inline hal_index_t MMapBottomSegment::getRightChildIndex(hal_size_t i) const
{
  assert(isLast() == false);
  MMapBottomSegment rightSeg(_genome, _index + 1);
  return rightSeg.getChildIndex(i);
}

}
#endif
