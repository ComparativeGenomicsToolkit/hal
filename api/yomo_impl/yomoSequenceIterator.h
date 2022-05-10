/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOSEQUENCEITERATOR_H
#define _YOMOSEQUENCEITERATOR_H

#include "halSequenceIterator.h"
#include "yomoExternalArray.h"
#include "yomoGenome.h"
#include "yomoSequence.h"
#include <H5Cpp.h>

namespace hal {

    class YomoSequenceIterator : public SequenceIterator {
      public:
        YomoSequenceIterator(YomoGenome *genome, hal_index_t index);
        ~YomoSequenceIterator();

        // SEQUENCE ITERATOR METHODS
        SequenceIteratorPtr clone() const;
        void toNext();
        void toPrev();
        bool atEnd() const;
        Sequence *getSequence();
        const Sequence *getSequence() const {
            return const_cast<YomoSequenceIterator *>(this)->getSequence();
        }
        bool equals(SequenceIteratorPtr other) const;

      private:
        YomoSequence _sequence;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
