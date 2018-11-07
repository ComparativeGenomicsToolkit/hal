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
    void toNext() { _sequence._data += 1; };
    void toPrev() { _sequence._data -= 1; };
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
