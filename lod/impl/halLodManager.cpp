/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halLodManager.h"
#include "halCLParser.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <limits>
#include <sstream>

#ifdef ENABLE_UDC
#include "udc2.h"
#endif

using namespace std;
using namespace hal;

// default to 5 days for now
const unsigned long LodManager::MaxAgeSec = 432000;

// specify upper limit of lods.
// (MUST MANUALLY KEEP CONSISTENT WITH global MaxLodToken variable in
// hal/lod/halLodInterpolate.py)
const string LodManager::MaxLodToken = "max";

LodManager::LodManager() : _options(NULL), _maxLodLowerBound((hal_size_t)numeric_limits<hal_index_t>::max()) {
    // FIXME: the way options work is weird.
}

LodManager::~LodManager() {
    for (AlignmentMap::iterator mapIt = _map.begin(); mapIt != _map.end(); ++mapIt) {
        if (mapIt->second.second.get() != NULL) {
            const_cast<Alignment *>(mapIt->second.second.get())->close();
        }
    }
}

void LodManager::loadLODFile(const string &lodPath, const CLParser *options) {
    _options = options;
    _map.clear();

#ifdef ENABLE_UDC
    char *cpath = const_cast<char *>(lodPath.c_str());

    size_t cbufSize = 0;
    char *cbuffer = udc2FileReadAll(cpath, NULL, 100000, &cbufSize);
    if (cbuffer == NULL) {
        throw hal_exception("Error udc-opening " + lodPath);
    }
    string cbufCpy(cbuffer);
    udc2FreeMem(cbuffer);
    stringstream ifile(cbufCpy);
#else
    ifstream ifile(lodPath.c_str());
#endif

    if (!ifile.good()) {
        stringstream emes;
        emes << "Error opening " << lodPath;
        throw hal_exception(emes.str());
    }

    string lineBuffer;
    hal_size_t minLen;
    string path;
    hal_size_t lineNum = 1;
    bool foundMax = false;
    _maxLodLowerBound = (hal_size_t)numeric_limits<hal_index_t>::max();
    while (ifile.good()) {
        getline(ifile, lineBuffer);
        stringstream ss(lineBuffer);
        ss >> minLen >> path;
        if (ifile.bad()) {
            throw hal_exception("Error parsing line " + std::to_string(lineNum) + " of " + lodPath);
        }
        if (lineBuffer.length() == 0) {
            continue;
        }
        string fullHalPath;
        if (foundMax == true) {
            throw hal_exception("Error on line " + std::to_string(lineNum) + " of " + lodPath + ": Limit token (" +
                                MaxLodToken + ") can only appear on" + " final line.");
        }
        if (path == MaxLodToken) {
            foundMax = true;
            _maxLodLowerBound = minLen;
            fullHalPath = MaxLodToken;
        } else {
            fullHalPath = resolvePath(lodPath, path);
        }
        _map.insert(pair<hal_size_t, PathAlign>(minLen, PathAlign(fullHalPath, AlignmentConstPtr())));
        ++lineNum;
    }

    checkMap(lodPath);
}

void LodManager::loadSingeHALFile(const string &halPath, const CLParser *options) {
    _options = options;
    _map.clear();
    _map.insert(pair<hal_size_t, PathAlign>(0, PathAlign(halPath, AlignmentConstPtr())));
    _maxLodLowerBound = (hal_size_t)numeric_limits<hal_index_t>::max();
    checkMap(halPath);
}

AlignmentConstPtr LodManager::getAlignment(hal_size_t queryLength, bool needDNA) {
    assert(_map.size() > 0);
    AlignmentMap::iterator mapIt;
    if (needDNA == true) {
        mapIt = _map.begin();
    } else {
        mapIt = _map.upper_bound(queryLength);
        --mapIt;
    }
    assert(mapIt->first <= queryLength);
    AlignmentConstPtr &alignment = mapIt->second.second;
    if (mapIt->first == _maxLodLowerBound) {
        throw hal_exception("Query length " + std::to_string(queryLength) + " above maximum LOD size of " +
                            std::to_string(getMaxQueryLength()));
    }
    if (alignment.get() == NULL) {
        alignment = AlignmentConstPtr(openHalAlignment(mapIt->second.first, _options));
        checkAlignment(mapIt->first, mapIt->second.first, alignment);
    }
    assert(mapIt->second.second.get() != NULL);
    return alignment;
}

bool LodManager::isLod0(hal_size_t queryLength) const {
    assert(_map.size() > 0);
    AlignmentMap::const_iterator mapIt = _map.upper_bound(queryLength);
    --mapIt;
    return mapIt == _map.begin();
}

string LodManager::resolvePath(const string &lodPath, const string &halPath) {
    assert(lodPath.empty() == false && halPath.empty() == false);
    if (halPath[0] == '/' || halPath.find(":/") != string::npos) {
        return halPath;
    }
    size_t sPos = lodPath.find_last_of('/');
    if (sPos == string::npos) {
        return halPath;
    }
    return lodPath.substr(0, sPos + 1) + halPath;
}

void LodManager::checkMap(const string &lodPath) {
    if (_map.size() == 0) {
        throw hal_exception("No entries were found in " + lodPath);
    }
    AlignmentMap::const_iterator mapIt = _map.begin();
    if (mapIt->first != 0) {
        throw hal_exception("No alignment with range value 0 found in " + lodPath + ". " +
                            "A record of the form \"0 pathToOriginalHALFile\" must be present");
    }
    if (_maxLodLowerBound == 0) {
        throw hal_exception("Maximum LOD query length must be > 0");
    }
}

void LodManager::checkAlignment(hal_size_t minQuery, const string &path, AlignmentConstPtr alignment) {
    if (alignment->getNumGenomes() == 0) {
        throw hal_exception("No genomes found in base alignment specified in " + path);
    }

#ifndef NDEBUG
    if (minQuery == 0) {
        vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
        string name = !leafNames.empty() ? leafNames[0] : alignment->getRootName();
        const Genome *genome = alignment->openGenome(name);

        bool seqFound = genome->containsDNAArray();
        alignment->closeGenome(genome);
        if (seqFound == false) {
            throw hal_exception("HAL file for highest level of detail (0) in genome " + name + " (" + path +
                                ") must contain DNA sequence information.");
        }
    }
#endif
}
