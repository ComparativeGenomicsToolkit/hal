/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoMetaData.h"
#include "halCommon.h"
#include "yomoCommon.h"
#include <cassert>
#include <iostream>

using namespace hal;
using namespace H5;
using namespace std;

YOMOMetaData::YOMOMetaData() : _parent(NULL) {
}

YOMOMetaData::YOMOMetaData(PortableH5Location *parent, const string &name) {
    open(parent, name);
}

YOMOMetaData::~YOMOMetaData() {
    write();
}

void YOMOMetaData::set(const string &key, const string &value) {
    if (has(key) == true) {
        _map[key] = value;
    } else {
        _map.insert(pair<string, string>(key, value));
    }
    _dirty = true;
}

const string &YOMOMetaData::get(const string &key) const {
    assert(has(key) == true);
    return _map.find(key)->second;
}

bool YOMOMetaData::has(const string &key) const {
    return _map.find(key) != _map.end();
}

const map<string, string> &YOMOMetaData::getMap() const {
    return _map;
}

// hack this in for compatibility for newer yomo which seems to have changed
// interface from object to location here.
#if H5_VERSION_GE(1, 8, 12) && H5_VERSION_LE(1, 10, 0)
#define ATTR_OP_PARAM__ H5Location
#else
#define ATTR_OP_PARAM__ H5Object
#endif
static void attr_operator(ATTR_OP_PARAM__ &loc /*in*/, const H5std_string attr_name /*in*/, void *operator_data /*in,out*/) {
    map<string, string> *attMap = static_cast<map<string, string> *>(operator_data);
    StrType vlsType(0, H5T_VARIABLE);
    Attribute attr = loc.openAttribute(attr_name);
    H5std_string strg;
    attr.read(vlsType, strg);
    attMap->insert(pair<string, string>(attr_name, strg));
}

void YOMOMetaData::open(PortableH5Location *parent, const string &name) {
    assert(parent != NULL);
    _map.clear();
    _parent = parent;
    _name = name;
    _dirty = false;

    try {
        YOMODisableExceptionPrinting prDisable;
        _group = parent->openGroup(name);
    } catch (Exception &e) {
        _group = parent->createGroup(name);
    }

    if (_group.getNumAttrs() > 0) {
        _group.iterateAttrs(attr_operator, NULL, (void *)&_map);
    }
    assert(_map.size() == (size_t)_group.getNumAttrs());
}

void YOMOMetaData::write() {
    if (!_dirty)
        return;

    for (map<string, string>::iterator i = _map.begin(); i != _map.end(); ++i) {
        StrType vlsType(0, H5T_VARIABLE);
        DataSpace attSpace(H5S_SCALAR);
        Attribute attr;
        // we delve into C api for this test
        if (H5Aexists(_group.getId(), i->first.c_str()) == true) {
            _group.removeAttr(i->first.c_str());
        }
        attr = _group.createAttribute(i->first, vlsType, attSpace);
        attr.write(vlsType, i->second);
    }
    _dirty = false;
}
