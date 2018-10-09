#ifndef _MMAPSEQUENCEITERATOR_H
#define _MMAPSEQUENCEITERATOR_H
#include "halSequenceIterator.h"
namespace hal {
class MMapSequenceIterator : public SequenceIterator {
public:
    // FIXME: these are all unimplemented currently, just want to get the compiler to stfu.
    MMapSequenceIterator(MMapGenome* genome, hal_index_t index) {};

    // SEQUENCE ITERATOR METHODS
    SequenceIteratorPtr copy() {return SequenceIteratorPtr(NULL);};
    SequenceIteratorConstPtr copy() const {return SequenceIteratorConstPtr(NULL);};
    void toNext() const {};
    void toPrev() const {};
    Sequence* getSequence() {return NULL;};
    const Sequence* getSequence() const {return NULL;};
    bool equals(SequenceIteratorConstPtr other) const {return false;};
};
}
#endif
