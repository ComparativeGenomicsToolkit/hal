/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halSegment.h"
#include "halGenome.h"
#include "halDnaIterator.h"

using namespace hal;

bool Segment::isMissingData(double nThreshold) const {
    DnaIteratorPtr dnaIt(getGenome()->getDnaIterator(getStartPosition()));
    size_t length = getLength();
    size_t maxNs = nThreshold * (double)length;
    size_t Ns = 0;
    char c;
    for (size_t i = 0; i < length; ++i, dnaIt->toRight()) {
        c = dnaIt->getBase();
        if (c == 'N' || c == 'n') {
            ++Ns;
        }
        if (Ns > maxNs) {
            return true;
        }
        if ((length - i) < (maxNs - Ns)) {
            return false;
        }
    }
    return false;
}

