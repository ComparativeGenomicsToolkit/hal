#ifndef _MMAPSEQUENCEITERATOR_H
#define _MMAPSEQUENCEITERATOR_H
#include "halSequenceIterator.h"
namespace hal {
// FIXME: the fact that _data is being moved but _index is not moved is confusing;
// why does this class have an _index?
class MMapSequenceIterator : public SequenceIterator {
public:
    MMapSequenceIterator(MMapGenome* genome, hal_index_t index) :
        _genome(genome),
        _sequence(_genome, _genome->getSequenceData(index)),
        _index(index)  {
    };

    // SEQUENCE ITERATOR METHODS
    SequenceIteratorPtr clone() const {
        return SequenceIteratorPtr(new MMapSequenceIterator(_genome, _index));
    }
    void toNext() {
        _index++;
        _sequence._data += 1;
    }
    void toPrev() {
        _index--;
        _sequence._data -= 1;
    }
    bool atEnd() const {
        return (_sequence._data < _genome->getSequenceData(0))
            or (_sequence._data >= _genome->getSequenceData(_genome->getNumSequences()));
    }
    const Sequence* getSequence() const {
        return &_sequence;
    }
    bool equals(SequenceIteratorPtr other) const {
        const MMapSequenceIterator* mmapOther = reinterpret_cast<
            const MMapSequenceIterator*>(other.get());
        assert(_sequence.getGenome() == mmapOther->_sequence.getGenome());
        return _sequence.getArrayIndex() == mmapOther->_sequence.getArrayIndex();
    };

private:
    MMapGenome *_genome;
    MMapSequence _sequence;
    hal_index_t _index;
};
}
#endif
