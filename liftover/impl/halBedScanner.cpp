/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "halBedScanner.h"

using namespace std;
using namespace hal;

BedScanner::BedScanner() : _bedStream(NULL) {
}

BedScanner::~BedScanner() {
}

void BedScanner::scan(const string &bedPath, int bedType) {
    assert(_bedStream == NULL);
    _bedStream = new ifstream(bedPath.c_str());
    try {
        scan(_bedStream, bedType);
    } catch (hal_exception &e) {
        delete _bedStream;
        _bedStream = NULL;
        throw hal_exception(string(e.what()) + " in file " + bedPath);
    }

    delete _bedStream;
    _bedStream = NULL;
}

void BedScanner::scan(istream *is, int bedType) {
    visitBegin();
    _bedStream = is;
    if (_bedStream->bad()) {
        throw hal_exception("Error reading bed input stream");
    }
    string lineBuffer;
    _lineNumber = 0;
    try {
        skipWhiteSpaces(_bedStream);
        while (_bedStream->good()) {
            ++_lineNumber;
            _bedLine.read(*_bedStream, lineBuffer, bedType);
            visitLine();
            skipWhiteSpaces(_bedStream);
        }
    } catch (hal_exception &e) {
        throw hal_exception(string(e.what()) + " in input bed line " + std::to_string(_lineNumber));
    }
    visitEOF();
    _bedStream = NULL;
}

size_t BedScanner::getNumColumns(const string &bedLine) {
    return chopString(bedLine, "\t").size();
}
void BedScanner::visitBegin() {
}

void BedScanner::visitLine() {
}

void BedScanner::visitEOF() {
}

void BedScanner::skipWhiteSpaces(istream *bedStream) {
    while (bedStream->good() && std::isspace((char)bedStream->peek())) {
        bedStream->get();
    }
}
