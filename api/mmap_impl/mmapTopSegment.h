#ifndef _MMAPTOPSEGMENT_H
#define _MMAPTOPSEGMENT_H
#include "halTopSegment.h"
#include "halGenome.h"
#include "mmapGenome.h"
#include "mmapTopSegmentData.h"
#include "cassert"

namespace hal {
class MMapTopSegment : public TopSegment
{
    public:
    MMapTopSegment(const MMapGenome *genome, hal_index_t arrayIndex) :
        _index(arrayIndex) {
            _genome = const_cast<MMapGenome *>(genome);
            _data = _genome->getTopSegmentPointer(arrayIndex);
        };

    // SEGMENT INTERFACE
    void setArrayIndex( Genome *genome, hal_index_t arrayIndex) {
        _genome = dynamic_cast<MMapGenome *>(genome);
        _data = _genome->getTopSegmentPointer(arrayIndex);
        _index = arrayIndex;
    }
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

    // TOP SEGMENT INTERFACE
    hal_index_t getParentIndex() const { return _data->getParentIndex(); };
    bool hasParent() const;
    void setParentIndex(hal_index_t parIdx) { _data->setParentIndex(parIdx); };
    bool getParentReversed() const { return _data->getReversed(); };
    void setParentReversed(bool isReversed) { _data->setReversed(isReversed); };
    hal_index_t getBottomParseIndex() const { return _data->getBottomParseIndex(); };
    void setBottomParseIndex(hal_index_t botParseIdx) { _data->setBottomParseIndex(botParseIdx); };
    hal_offset_t getBottomParseOffset() const;
    bool hasParseDown() const;
    hal_index_t getNextParalogyIndex() const { return _data->getNextParalogyIndex(); }
    bool hasNextParalogy() const;
    void setNextParalogyIndex(hal_index_t parIdx) { _data->setNextParalogyIndex(parIdx); };
    hal_index_t getLeftParentIndex() const;
    hal_index_t getRightParentIndex() const;
    bool isCanonicalParalog() const;

    private:
    MMapTopSegmentData *_data;
    MMapGenome *_genome;
    hal_index_t _index;
};

inline const Genome* MMapTopSegment::getGenome() const
{
  return _genome;
}

inline Genome* MMapTopSegment::getGenome()
{
  return _genome;
}

inline hal_index_t MMapTopSegment::getEndPosition() const
{
  return getStartPosition() + (hal_index_t)(getLength() - 1);
}

inline hal_size_t MMapTopSegment::getLength() const
{
  return (_data + 1)->getStartPosition() - _data->getStartPosition();
}

inline const Sequence* MMapTopSegment::getSequence() const
{
  return _genome->getSequenceBySite(getStartPosition());
}

inline bool MMapTopSegment::hasParseDown() const
{
  return getBottomParseIndex() != NULL_INDEX;
}

inline bool MMapTopSegment::hasNextParalogy() const
{
  return getNextParalogyIndex() != NULL_INDEX;
}

inline bool MMapTopSegment::hasParent() const
{
  return getParentIndex() != NULL_INDEX;
}

inline hal_index_t MMapTopSegment::getArrayIndex() const
{
  return _index;
}

inline bool MMapTopSegment::leftOf(hal_index_t genomePos) const
{
  return getEndPosition() < genomePos;
}

inline bool MMapTopSegment::rightOf(hal_index_t genomePos) const
{
  return getStartPosition() > genomePos;
}

inline bool MMapTopSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

inline bool MMapTopSegment::isFirst() const
{
  assert(getSequence() != NULL);
  return _index == 0 || 
     _index == (hal_index_t)getSequence()->getTopSegmentArrayIndex();
}

inline bool MMapTopSegment::isLast() const
{
  assert(getSequence() != NULL);
  return _index == _genome->getNumTopSegments() || 
     _index == getSequence()->getTopSegmentArrayIndex() +
     (hal_index_t)getSequence()->getNumTopSegments() - 1;
}

inline bool MMapTopSegment::isTop() const
{
  return true;
}

inline hal_size_t MMapTopSegment::getMappedSegments(
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

inline hal_index_t MMapTopSegment::getLeftParentIndex() const
{
  assert(isFirst() == false);
  MMapTopSegment leftSeg(_genome, _index - 1);
  return leftSeg.getParentIndex();
}

inline hal_index_t MMapTopSegment::getRightParentIndex() const
{
  assert(isLast() == false);
  MMapTopSegment rightSeg(_genome, _index + 1);
  return rightSeg.getParentIndex();
}
}
#endif
