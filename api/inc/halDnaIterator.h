/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDNAITERATOR_H
#define _HALDNAITERATOR_H
#include "halDefs.h"
#include "halDnaDriver.h"
#include "halGenome.h"

namespace hal {

    /**
     * Interface for general dna iterator
     */
    class DnaIterator {
      public:
        DnaIterator(Genome *genome, DnaAccessPtr dnaAccess, hal_index_t index)
            : _genome(genome), _dnaAccess(dnaAccess), _index(index), _reversed(false) {
        }

        /** Destructor */
        virtual ~DnaIterator() {
        }

        /** Ensure cache write data is flushed if needed.  This should be called
         * after a writing loop or an error will be generated.  This is not done
         * automatically on destruct, as error will be caught in DnaAccess.
         */
        void flush() {
            _dnaAccess->flush();
        }

        /** Get the DNA character at this position (if revsersed is set
        * the reverse compelement is returned */
        char getBase() const;

        /** Set the DNA character at this position (if revsersed is set
         * the reverse compelement is stored
         * @param c DNA character to set */
        void setBase(char c);

        /** Move to previous position (equiv. to toRight if reversed)*/
        void toLeft() {
            _reversed ? ++_index : --_index;
        }

        /** Move to next position (equiv. to toLeft if reversed)*/
        void toRight() {
            _reversed ? --_index : ++_index;
        }

        /** Jump to any point on the genome (can lead to
         * inefficient paging from disk if used irresponsibly)
         * @param index position in array to jump to */
        void jumpTo(hal_size_t index) {
            _index = index;
        }

        /** has the iterator reach the end of the traversal in the direction of
         * movement? */
        bool atEnd() const {
            if (not _reversed) {
                return _index >= (hal_index_t)_genome->getSequenceLength();
            } else {
                return _index < 0;
            }
        }

        /** switch to base's reverse complement */
        void toReverse() {
            _reversed = !_reversed;
        }

        /** Check whether iterator is on base's complement */
        bool getReversed() const {
            return _reversed;
        }

        /** Set the iterator's reverse complement status */
        void setReversed(bool reversed) {
            _reversed = reversed;
        }

        /** Get the containing (read-only) genome */
        const Genome *getGenome() const {
            return _genome;
        }

        /** Get the containing genome */
        Genome *getGenome() {
            return _genome;
        }

        /** Get the containing sequence */
        const Sequence *getSequence() const {
            return _genome->getSequenceBySite(_index);
        }

        /** Get the index of the base in the dna array */
        hal_index_t getArrayIndex() const {
            return _index;
        }

        /* read a DNA string */
        void readString(std::string &outString, hal_size_t length);

        /* write a DNA string */
        void writeString(const std::string &inString, hal_size_t length);

        /** Compare (array indexes) of two iterators */
        bool equals(DnaIteratorPtr &other) const;

        /** Compare (array indexes) of two iterators */
        bool leftOf(DnaIteratorPtr &other) const;

      private:
        bool inRange() const {
            return (_index >= 0) && (_index < (hal_index_t)_genome->getSequenceLength());
        }

        Genome *_genome;
        DnaAccessPtr _dnaAccess;
        hal_index_t _index;
        bool _reversed;
    };

    inline char DnaIterator::getBase() const {
        assert(inRange());
        char c = _dnaAccess->getBase(_index);
        if (_reversed) {
            c = reverseComplement(c);
        }
        return c;
    }

    inline void DnaIterator::setBase(char c) {
        if (not inRange()) {
            throw hal_exception("Trying to set character out of range");
        } else if (not isNucleotide(c)) {
            throw hal_exception(std::string("Trying to set invalid character: ") + c);
        }
        if (_reversed) {
            c = reverseComplement(c);
        }
        _dnaAccess->setBase(_index, c);
        assert(getBase() == !_reversed ? c : reverseComplement(c));
    }

    inline bool DnaIterator::equals(DnaIteratorPtr &other) const {
        const DnaIterator *mmOther = reinterpret_cast<const DnaIterator *>(other.get());
        assert(_genome == mmOther->_genome);
        return _index == mmOther->_index;
    }

    inline bool DnaIterator::leftOf(DnaIteratorPtr &other) const {
        const DnaIterator *mmOther = reinterpret_cast<const DnaIterator *>(other.get());
        assert(_genome == mmOther->_genome);
        return _index < mmOther->_index;
    }

    inline void DnaIterator::readString(std::string &outString, hal_size_t length) {
        assert(length == 0 || inRange() == true);
        outString.resize(length);

        for (hal_size_t i = 0; i < length; ++i) {
            outString[i] = getBase();
            toRight();
        }
    }

    inline void DnaIterator::writeString(const std::string &inString, hal_size_t length) {
        assert(length == 0 || inRange());
        for (hal_size_t i = 0; i < length; ++i) {
            setBase(inString[i]);
            toRight();
        }
        flush();
    }
}

#endif
// Local Variables:
// mode: c++
// End:
