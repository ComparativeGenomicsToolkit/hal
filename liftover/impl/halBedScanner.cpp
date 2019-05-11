/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
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

void BedScanner::scan(const string &bedPath, int bedVersion, const locale *inLocale) {
    assert(_bedStream == NULL);
    _bedStream = new ifstream(bedPath.c_str());
    if (inLocale != NULL) {
        _bedStream->imbue(*inLocale);
    }

    try {
        scan(_bedStream, bedVersion);
    } catch (hal_exception &e) {
        delete _bedStream;
        _bedStream = NULL;
        throw hal_exception(string(e.what()) + " in file " + bedPath);
    }

    delete _bedStream;
    _bedStream = NULL;
}

void BedScanner::scan(istream *is, int bedVersion, const locale *inLocale) {
    visitBegin();
    _bedStream = is;
    _bedVersion = bedVersion;
    if (inLocale != NULL) {
        _bedStream->imbue(*inLocale);
    }

    if (_bedVersion == -1) {
        _bedVersion = getBedVersion(is);
    }

    if (_bedStream->bad()) {
        throw hal_exception("Error reading bed input stream");
    }
    string lineBuffer;
    _lineNumber = 0;
    try {
        skipWhiteSpaces(_bedStream, inLocale);
        while (_bedStream->good()) {
            ++_lineNumber;
            _bedLine.read(*_bedStream, _bedVersion, lineBuffer);
            visitLine();
            skipWhiteSpaces(_bedStream, inLocale);
        }
    } catch (hal_exception &e) {
        throw hal_exception(string(e.what()) + " in input bed line " + std::to_string(_lineNumber));
    }
    visitEOF();
    _bedStream = NULL;
}

int BedScanner::getBedVersion(istream *bedStream, const locale *inLocale) {
    assert(bedStream != &cin);
    if (bedStream->bad()) {
        throw hal_exception("Error reading bed input stream");
    }
    if (inLocale != NULL) {
        bedStream->imbue(*inLocale);
    }

    string lineBuffer;
    BedLine bedLine;
    int version = 12;
    streampos pos = bedStream->tellg();
    for (; version > 3; --version) {
        try {
            bedStream->clear();
            bedStream->seekg(pos);
            skipWhiteSpaces(bedStream, inLocale);
            *bedStream >> std::skipws;
            bedLine.read(*bedStream, version, lineBuffer);
            break;
        } catch (...) {
            if (version == 12) {
                version = 10;
            }
        }
    }
    bedStream->clear();
    bedStream->seekg(pos);
    assert(!bedStream->bad());
    return version;
}

size_t BedScanner::getNumColumns(const string &bedLine, const locale *inLocale) {
    stringstream ss(bedLine);
    if (inLocale != NULL) {
        ss.imbue(*inLocale);
    }
    size_t c = 0;
    string buffer;
    while (ss.good()) {
        buffer.clear();
        ss >> buffer;
        if (!buffer.empty()) {
            ++c;
        }
    }
    return c;
}
void BedScanner::visitBegin() {
}

void BedScanner::visitLine() {
}

void BedScanner::visitEOF() {
}

void BedScanner::skipWhiteSpaces(istream *bedStream, const locale *inLocale) {
    locale defaultLocale;
    const locale &myLocale = inLocale == NULL ? defaultLocale : *inLocale;
    while (bedStream->good() && std::isspace((char)bedStream->peek(), myLocale)) {
        bedStream->get();
    }
}
