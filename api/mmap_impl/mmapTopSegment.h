#ifndef _MMAPTOPSEGMENT_H
#define _MMAPTOPSEGMENT_H
#include "cassert"
#include "halGenome.h"
#include "halTopSegment.h"
#include "mmapGenome.h"
#include "mmapTopSegmentData.h"

namespace hal {
    class MMapTopSegment : public TopSegment {
      public:
        MMapTopSegment(MMapGenome *genome, hal_index_t arrayIndex)
            : TopSegment(genome, arrayIndex), _data(getMMapGenome()->getTopSegmentPointer(arrayIndex)) {
        }

        // SEGMENT INTERFACE
        void setArrayIndex(Genome *genome, hal_index_t arrayIndex) {
            _genome = genome;
            _data = getMMapGenome()->getTopSegmentPointer(arrayIndex);
            _index = arrayIndex;
        }
        const Sequence *getSequence() const;
        hal_index_t getStartPosition() const {
            return _data->getStartPosition();
        };
        hal_index_t getEndPosition() const;
        hal_size_t getLength() const;
        void setCoordinates(hal_index_t startPos, hal_size_t length);
        hal_index_t getArrayIndex() const;
        bool isFirst() const;
        bool isLast() const;
        std::ostream &print(std::ostream &os) const;

        // TOP SEGMENT INTERFACE
        hal_index_t getParentIndex() const {
            return _data->getParentIndex();
        };
        bool hasParent() const;
        void setParentIndex(hal_index_t parIdx) {
            _data->setParentIndex(parIdx);
        };
        bool getParentReversed() const {
            return _data->getReversed();
        };
        void setParentReversed(bool isReversed) {
            _data->setReversed(isReversed);
        };
        hal_index_t getBottomParseIndex() const {
            return _data->getBottomParseIndex();
        };
        void setBottomParseIndex(hal_index_t botParseIdx) {
            _data->setBottomParseIndex(botParseIdx);
        };
        hal_offset_t getBottomParseOffset() const;
        bool hasParseDown() const;
        hal_index_t getNextParalogyIndex() const {
            return _data->getNextParalogyIndex();
        }
        bool hasNextParalogy() const;
        void setNextParalogyIndex(hal_index_t parIdx) {
            _data->setNextParalogyIndex(parIdx);
        };
        hal_index_t getLeftParentIndex() const;
        hal_index_t getRightParentIndex() const;
        bool isCanonicalParalog() const;

      private:
        MMapGenome *getMMapGenome() const {
            return static_cast<MMapGenome *>(_genome);
        }
        MMapTopSegmentData *_data;
    };

    inline hal_index_t MMapTopSegment::getEndPosition() const {
        return getStartPosition() + (hal_index_t)(getLength() - 1);
    }

    inline hal_size_t MMapTopSegment::getLength() const {
        return (_data + 1)->getStartPosition() - _data->getStartPosition();
    }

    inline const Sequence *MMapTopSegment::getSequence() const {
        return _genome->getSequenceBySite(getStartPosition());
    }

    inline bool MMapTopSegment::hasParseDown() const {
        return getBottomParseIndex() != NULL_INDEX;
    }

    inline bool MMapTopSegment::hasNextParalogy() const {
        return getNextParalogyIndex() != NULL_INDEX;
    }

    inline bool MMapTopSegment::hasParent() const {
        return getParentIndex() != NULL_INDEX;
    }

    inline hal_index_t MMapTopSegment::getArrayIndex() const {
        return _index;
    }

    inline bool MMapTopSegment::isFirst() const {
        assert(getSequence() != NULL);
        return _index == 0 || _index == (hal_index_t)getSequence()->getTopSegmentArrayIndex();
    }

    inline bool MMapTopSegment::isLast() const {
        assert(getSequence() != NULL);
        return _index == _genome->getNumTopSegments() ||
               _index == getSequence()->getTopSegmentArrayIndex() + (hal_index_t)getSequence()->getNumTopSegments() - 1;
    }

    inline hal_index_t MMapTopSegment::getLeftParentIndex() const {
        assert(isFirst() == false);
        MMapTopSegment leftSeg(getMMapGenome(), _index - 1);
        return leftSeg.getParentIndex();
    }

    inline hal_index_t MMapTopSegment::getRightParentIndex() const {
        assert(isLast() == false);
        MMapTopSegment rightSeg(getMMapGenome(), _index + 1);
        return rightSeg.getParentIndex();
    }
}
#endif
// Local Variables:
// mode: c++
// End:
