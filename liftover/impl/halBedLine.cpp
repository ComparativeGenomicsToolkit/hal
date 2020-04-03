/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "halBedLine.h"
#include "halCommon.h"

using namespace std;
using namespace hal;

BedLine::BedLine() : _start(NULL_INDEX), _end(NULL_INDEX), _strand('+'), _bedType(NULL_INDEX), _srcStart(NULL_INDEX) {
}

BedLine::~BedLine() {
}

/* bedType is zero or the number of standard bed columns.  All others
 * are saved as extra */
istream &BedLine::read(istream &is, string &lineBuffer, int bedType) {
    _bedType = bedType;
    std::getline(is, lineBuffer);
    std::vector<std::string> row = chopString(lineBuffer, "\t");
    if (row.size() < 3) {
        throw hal_exception("Expected at least three columns in BED record: " + lineBuffer);
    }
    if (_bedType == 0) {
        _bedType = min(int(row.size()), 12);
    }
    _chrName = row[0];
    _start = strToInt(row[1]);
    _end = strToInt(row[2]);
    if (_start >= _end) {
        throw hal_exception("Error zero or negative length BED range: " + lineBuffer);
    }
    if (_bedType > 3) {
        _name = row[3];
    }
    if (_bedType > 4) {
        _score = strToInt(row[4]);
    }
    if (_bedType > 5) {
        _strand = row[5][0];
        if (_strand != '.' && _strand != '+' && _strand != '-') {
            throw hal_exception("Strand character must be + or - or ." + lineBuffer);
        }
    }
    if (_bedType > 6) {
        _thickStart = strToInt(row[6]);
    }
    if (_bedType > 7) {
        _thickEnd = strToInt(row[7]);
    }
   if (_bedType > 8) {
        vector<string> rgbTokens = chopString(row[8], ",");
        if (rgbTokens.size() > 3 || rgbTokens.size() == 0) {
            throw hal_exception("Error parsing BED itemRGB: " + lineBuffer);
        }
        _itemR = strToInt(rgbTokens[0]);
        _itemG = _itemB = _itemR;
        if (rgbTokens.size() > 1) {
            _itemG = strToInt(rgbTokens[1]);
        }
        if (rgbTokens.size() == 3) {
            _itemB= strToInt(rgbTokens[2]);;
        }
    }
    if (_bedType > 9) {
        if (_bedType < 12) {
            throw hal_exception("Error parsing BED, insufficient columns for blocks: " + lineBuffer);
        }
        size_t numBlocks = strToInt(row[9]);
        vector<string> blockSizes = chopString(row[10], ",");
        if (blockSizes.size() != numBlocks) {
            throw hal_exception("Error parsing BED blockSizes: " + lineBuffer);
        }
        vector<string> blockStarts = chopString(row[11], ",");
        if (blockStarts.size() != numBlocks) {
            throw hal_exception("Error parsing BED blockStarts: " + lineBuffer);
        }
        _blocks.resize(numBlocks);
        for (size_t i = 0; i < numBlocks; ++i) {
            _blocks[i]._length = strToInt(blockSizes[i]);
            _blocks[i]._start = strToInt(blockStarts[i]);
            if (_start + _blocks[i]._start + _blocks[i]._length > _end) {
                throw hal_exception("Error BED block out of range: " + lineBuffer);
            }
        }
    }
    for (int i = _bedType; i < row.size(); i++) {
        _extra.push_back(row[i]);
    }
    return is;
}

ostream &BedLine::write(ostream &os) {
    os << _chrName << '\t' << _start << '\t' << _end;

    if (_bedType > 3) {
        os << '\t' << _name;
    }
    if (_bedType > 4) {
        os << '\t' << _score;
    }
    if (_bedType > 5) {
        os << '\t' << _strand;
    }
    if (_bedType > 6) {
        os << '\t' << _thickStart;
    }
    if (_bedType > 7) {
        os << '\t' << _thickEnd;
    }
    if (_bedType > 8) {
        os << '\t' << _itemR << ',' << _itemG << ',' << _itemB;
    }
    if (_bedType > 9) {
        os << '\t' << _blocks.size();
        for (size_t i = 0; i < _blocks.size(); ++i) {
            if (i == 0) {
                os << '\t';
            } else {
                os << ',';
            }
            os << _blocks[i]._length;
        }
        for (size_t i = 0; i < _blocks.size(); ++i) {
            if (i == 0) {
                os << '\t';
            } else {
                os << ',';
            }
            os << _blocks[i]._start;
        }
    }

    for (size_t i = 0; i < _extra.size(); ++i) {
        os << '\t' << _extra[i];
    }

    os << '\n';
    return os;
}

void BedLine::expandToBed12() {
    if (_bedType <= 3) {
        _name = "";
    }
    if (_bedType <= 4) {
        _score = 0;
    }
    if (_bedType <= 5) {
        _strand = '+';
    }
    if (_bedType <= 6) {
        _thickStart = _start;
    }
    if (_bedType <= 7) {
        _thickEnd = _end;
    }
    if (_bedType <= 8) {
        _itemR = _itemG = _itemB = 0;
    }
    if (_bedType <= 9) {
        _blocks.resize(1);
        _blocks[0]._start = 0;
        _blocks[0]._length = _end - _start;
    }
    _bedType = 12;
}

bool BedBlock::operator<(const BedBlock &other) const {
    bool isLess = _start < other._start;
    return isLess;
}

bool BedLineLess::operator()(const BedLine &b1, const BedLine &b2) const {
    if (b1._chrName < b2._chrName) {
        return true;
    } else if (b1._chrName == b2._chrName) {
        if (b1._start < b2._start) {
            return true;
        } else if (b1._start == b2._start) {
            if (b1._end < b2._end) {
                return true;
            } else if (b1._end == b2._end) {
                return b1._strand == '+' && b2._strand != '+';
            }
        }
    }
    return false;
}

bool BedLineSrcLess::operator()(const BedLine &b1, const BedLine &b2) const {
    return b1._srcStart < b2._srcStart;
}

ostream &BedLine::writePSL(ostream &os, bool prefixWithName) {
    assert(_psl.size() == 1);
    const PSLInfo &psl = _psl[0];
    assert(_blocks.size() == psl._qBlockStarts.size());
    assert(_blocks.size() > 0);
    assert(_srcStart >= (hal_index_t)psl._qChromOffset);
    if (validatePSL() == false) {
        throw hal_exception("Internal error: PSL does not validate");
    }

    if (prefixWithName == true) {
        os << _name << '\t';
    }
    os << psl._matches << '\t' << psl._misMatches << '\t' << psl._repMatches << '\t' << psl._nCount << '\t' << psl._qNumInsert
       << '\t' << psl._qBaseInsert << '\t' << psl._tNumInsert << '\t' << psl._tBaseInsert << '\t' << psl._qStrand << _strand
       << '\t' << psl._qSeqName << '\t' << psl._qSeqSize << '\t' << (_srcStart - psl._qChromOffset) << '\t'
       << (psl._qEnd - psl._qChromOffset) << '\t' << _chrName << '\t' << psl._tSeqSize << '\t' << _start << '\t' << _end << '\t'
       << _blocks.size() << '\t';

    for (size_t i = 0; i < _blocks.size(); ++i) {
        os << _blocks[i]._length << ',';
    }
    os << '\t';

    for (size_t i = 0; i < psl._qBlockStarts.size(); ++i) {
        assert(psl._qBlockStarts[i] >= (hal_index_t)psl._qChromOffset);
        hal_index_t start = psl._qBlockStarts[i] - psl._qChromOffset;
        if (psl._qStrand == '-') {
            start = psl._qSeqSize - start - _blocks[i]._length;
        }
        os << start << ',';
    }
    os << '\t';

    for (size_t i = 0; i < _blocks.size(); ++i) {
        hal_index_t start = _blocks[i]._start + _start;
        if (_strand == '-') {
            start = psl._tSeqSize - start - _blocks[i]._length;
        }
        os << start << ',';
    }

    os << '\n';
    return os;
}

bool BedLine::validatePSL() const {
    if (_psl.size() != 1) {
        assert(false);
        return false;
    }
    const PSLInfo &psl = _psl[0];

    if (_blocks.size() < 1) {
        assert(false);
        return false;
    }
    if (_blocks.size() != psl._qBlockStarts.size()) {
        assert(false);
        return false;
    }

    if (psl._qBlockStarts.size() != _blocks.size()) {
        assert(false);
        return false;
    }

    hal_size_t totBlockLen = 0;
    for (size_t i = 0; i < _blocks.size(); ++i) {
        totBlockLen += _blocks[i]._length;
    }

    if (totBlockLen != psl._matches + psl._misMatches + psl._repMatches + psl._nCount) {
        assert(false);
        return false;
    }

    if (totBlockLen + psl._qBaseInsert != psl._qEnd - _srcStart) {
        assert(false);
        return false;
    }

    if (totBlockLen + psl._tBaseInsert != (hal_size_t)_end - _start) {
        assert(false);
        return false;
    }

    if (_strand != '-') {
        if (_blocks[0]._start != 0) {
            assert(false);
            return false;
        }
        if (_blocks.back()._start + _blocks.back()._length + _start != _end) {
            assert(false);
            return false;
        }
    } else {
        if (_blocks.back()._start != 0) {
            assert(false);
            return false;
        }
        if (_blocks[0]._start + _blocks[0]._length + _start != _end) {
            assert(false);
            return false;
        }
    }

    if (psl._qStrand != '-') {
        if (psl._qBlockStarts[0] - psl._qChromOffset != _srcStart - psl._qChromOffset) {
            assert(false);
            return false;
        }
        if (psl._qBlockStarts.back() + _blocks.back()._length != psl._qEnd) {
            assert(false);
            return false;
        }
    } else {
        if (psl._qBlockStarts.back() - psl._qChromOffset != _srcStart - psl._qChromOffset) {
            assert(false);
            return false;
        }
        if (psl._qBlockStarts[0] + _blocks[0]._length != psl._qEnd) {
            assert(false);
            return false;
        }
    }

    return true;
}
