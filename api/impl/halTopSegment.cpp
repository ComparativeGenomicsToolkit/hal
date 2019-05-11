#include "halTopSegment.h"
#include "halDnaIterator.h"

using namespace hal;

void TopSegment::getString(std::string &outString) const {
    DnaIteratorPtr dnaIt(getGenome()->getDnaIterator(getStartPosition()));
    dnaIt->readString(outString, getLength());
}
