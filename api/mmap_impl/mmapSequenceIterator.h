#ifndef _MMAPSEQUENCEITERATOR_H
#define _MMAPSEQUENCEITERATOR_H
#include "halSequenceIterator.h"
namespace hal {
class MMapSequenceIterator : public SequenceIterator {
public:
    MMapSequenceIterator(MMapGenome* genome, hal_index_t index) :
        _genome(genome),
        _sequence(_genome, _genome->getSequenceData(index)),
        _index(index)  {
    };

    // SEQUENCE ITERATOR METHODS
    SequenceIteratorPtr clone() const { return SequenceIteratorPtr(new MMapSequenceIterator(_genome, _index)); };
    void toNext() { _index++; _sequence._data = _genome->getSequenceData(_index); };
    void toPrev() { _index--; _sequence._data = _genome->getSequenceData(_index); };
    const Sequence* getSequence() const { return &_sequence; };
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
