/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halWiggleLoader.h"
#include "halWiggleLiftover.h"
#include <cassert>
#include <deque>

using namespace std;
using namespace hal;

WiggleLoader::WiggleLoader() {
}

WiggleLoader::~WiggleLoader() {
}

void WiggleLoader::load(const Alignment *alignment, const Genome *genome, istream *inputFile, WiggleTiles<double> *vals) {
    _alignment = AlignmentConstPtr(alignment);
    _srcGenome = genome;
    _srcSequence = NULL;
    _vals = vals;
    scan(inputFile);
}

void WiggleLoader::visitHeader() {
    _srcSequence = _srcGenome->getSequence(_sequenceName);
    if (_srcSequence == NULL) {
        throw hal_exception("Sequence " + _sequenceName + " not found in genome " + _srcGenome->getName());
    }
}

void WiggleLoader::visitLine() {
    if (_srcSequence == NULL) {
        throw hal_exception("Missing Wig header");
    }

    hal_index_t absFirst = _first + _srcSequence->getStartPosition();
    hal_index_t absLast = _last + _srcSequence->getStartPosition();

    for (hal_index_t absPos = absFirst; absPos <= absLast; ++absPos) {
        _vals->set(absPos, _value);
    }
}
