/* Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoGenome.h"
#include "H5Cpp.h"
#include "halBottomSegmentIterator.h"
#include "halColumnIterator.h"
#include "halDnaIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halRearrangement.h"
#include "halTopSegmentIterator.h"
#include "yomoBottomSegment.h"
#include "yomoCommon.h"
#include "yomoDnaDriver.h"
#include "yomoSequence.h"
#include "yomoSequenceIterator.h"
#include "yomoTopSegment.h"
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace hal;
using namespace std;
using namespace H5;

const string YomoGenome::dnaArrayName = "DNA_ARRAY";
const string YomoGenome::topArrayName = "TOP_ARRAY";
const string YomoGenome::bottomArrayName = "BOTTOM_ARRAY";
const string YomoGenome::sequenceIdxArrayName = "SEQIDX_ARRAY";
const string YomoGenome::sequenceNameArrayName = "SEQNAME_ARRAY";
const string YomoGenome::metaGroupName = "Meta";
const string YomoGenome::rupGroupName = "Rup";
const double YomoGenome::dnaChunkScale = 10.;

YomoGenome::YomoGenome(const string &name, YomoAlignment *alignment, PortableH5Location *h5Parent,
                       const DSetCreatPropList &dcProps, bool inMemory)
    : Genome(alignment, name), _alignment(alignment), _h5Parent(h5Parent), _name(name), _numChildrenInBottomArray(0),
      _totalSequenceLength(0), _numChunksInArrayBuffer(inMemory ? 0 : 1) {
    _dcprops.copy(dcProps);
    assert(!name.empty());
    assert(alignment != NULL && h5Parent != NULL);

    try {
        YOMODisableExceptionPrinting prDisable;
        _group = h5Parent->openGroup(name);
    } catch (Exception &e) {
        _group = h5Parent->createGroup(name);
    }
    read();
    _metaData = new YOMOMetaData(&_group, metaGroupName);
    _rup = new YOMOMetaData(&_group, rupGroupName);

    _totalSequenceLength = _dnaArray.getSize() * 2;
    if (_totalSequenceLength > 0 && _rup->get(rupGroupName) == "1") {
        _totalSequenceLength -= 1;
    } else if (_totalSequenceLength == 0 && _sequenceIdxArray.getSize() > 0) {
        YomoSequence lastSeq(this, &_sequenceIdxArray, &_sequenceNameArray, _sequenceNameArray.getSize() - 1);
        _totalSequenceLength = lastSeq.getEndPosition() + 1;
    }
}

YomoGenome::~YomoGenome() {
    delete _metaData;
    delete _rup;
    deleteSequenceCache();
}

// GENOME INTERFACE

void YomoGenome::setDimensions(const vector<Sequence::Info> &sequenceDimensions, bool storeDNAArrays) {
    _totalSequenceLength = 0;
    hal_size_t totalSeq = sequenceDimensions.size();
    hal_size_t maxName = 0;

    // Copy segment dimensions to use the external interface
    vector<Sequence::UpdateInfo> topDimensions;
    topDimensions.reserve(sequenceDimensions.size());
    vector<Sequence::UpdateInfo> bottomDimensions;
    bottomDimensions.reserve(sequenceDimensions.size());

    // Compute summary info from the list of sequence Dimensions
    for (vector<Sequence::Info>::const_iterator i = sequenceDimensions.begin(); i != sequenceDimensions.end(); ++i) {
        _totalSequenceLength += i->_length;
        maxName = max(static_cast<hal_size_t>(i->_name.length()), maxName);
        topDimensions.push_back(Sequence::UpdateInfo(i->_name, i->_numTopSegments));
        bottomDimensions.push_back(Sequence::UpdateInfo(i->_name, i->_numBottomSegments));
    }

    // Unlink the DNA and segment arrays if they exist (using
    // exceptions is the only way I know how right now).  Note that
    // the file needs to be refactored to take advantage of the new
    // space. (FIXME)
    try {
        YOMODisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(dnaArrayName);
        _group.unlink(dnaArrayName);
    } catch (H5::Exception &) {
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(sequenceIdxArrayName);
        _group.unlink(sequenceIdxArrayName);
    } catch (H5::Exception &) {
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(sequenceNameArrayName);
        _group.unlink(sequenceNameArrayName);
    } catch (H5::Exception &) {
    }

    if (_totalSequenceLength > 0 && storeDNAArrays) {
        hal_size_t arrayLength = _totalSequenceLength / 2;
        if (_totalSequenceLength % 2) {
            ++arrayLength;
            _rup->set(rupGroupName, "1");
        } else {
            _rup->set(rupGroupName, "0");
        }
        hsize_t chunk;
        _dcprops.getChunk(1, &chunk);
        // enalarge chunk size because dna bases are so much smaller
        // than segments.  (about 30x). we default to 10x enlargement
        // since the seem to compress about 3x worse.
        chunk *= dnaChunkScale;
        DSetCreatPropList dnaDC;
        dnaDC.copy(_dcprops);
        dnaDC.setChunk(1, &chunk);
        _dnaArray.create(&_group, dnaArrayName, dnaDataType(), arrayLength, &dnaDC, _numChunksInArrayBuffer);
        _dnaAccess = DnaAccessPtr(new YOMODnaAccess(this, &_dnaArray, 0));
    }
    if (totalSeq > 0) {
        _sequenceIdxArray.create(&_group, sequenceIdxArrayName, YomoSequence::idxDataType(), totalSeq + 1, &_dcprops,
                                 _numChunksInArrayBuffer);

        _sequenceNameArray.create(&_group, sequenceNameArrayName, YomoSequence::nameDataType(maxName + 1), totalSeq, &_dcprops,
                                  _numChunksInArrayBuffer);

        writeSequences(sequenceDimensions);
    }

    // Do the same as above for the segments.
    setGenomeTopDimensions(topDimensions);
    setGenomeBottomDimensions(bottomDimensions);

    reload();
}

void YomoGenome::updateTopDimensions(const vector<Sequence::UpdateInfo> &topDimensions) {
    loadSequencePosCache();
    loadSequenceNameCache();
    vector<Sequence::UpdateInfo>::const_iterator i;
    map<string, YomoSequence *>::iterator cacheIt;
    map<string, const Sequence::UpdateInfo *> inputMap;
    map<string, hal_size_t> currentTopD;
    // copy input into map, checking everything is already present
    for (i = topDimensions.begin(); i != topDimensions.end(); ++i) {
        const string &name = i->_name;
        cacheIt = _sequenceNameCache.find(name);
        if (cacheIt == _sequenceNameCache.end()) {
            throw hal_exception(string("Cannot update sequence ") + name + " because it is not present in "
                                                                           " genome " +
                                getName());
        }
        inputMap.insert(pair<string, const Sequence::UpdateInfo *>(name, &*i));
    }
    // keep a record of the number of segments in each existing
    // segment (these can get muddled as we add the new ones in the next
    // loop to be sure by getting them in one shot)
    // Note to self: iterating the map in this way skips zero-length
    // sequences (which are in the separate vector).  This is fine
    // here since we will never update them, but seems like it could be
    // dangerous if something were to change
    map<hal_size_t, YomoSequence *>::iterator posCacheIt;
    map<string, const Sequence::UpdateInfo *>::iterator inputIt;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        YomoSequence *sequence = posCacheIt->second;
        inputIt = inputMap.find(sequence->getName());
        if (inputIt == inputMap.end()) {
            currentTopD.insert(pair<string, hal_size_t>(sequence->getName(), sequence->getNumTopSegments()));
        }
    }
    // scan through existing sequences, updating as necessary
    // build summary of all new and unchanged dimensions in newDimensions
    // Note to self: iterating the map in this way skips zero-length
    // sequences (which are in the separate vector).  This is fine
    // here since we will never update them, but seems like it could be
    // dangerous if something were to change
    map<string, hal_size_t>::iterator currentIt;
    vector<Sequence::UpdateInfo> newDimensions;
    Sequence::UpdateInfo newInfo;
    hal_size_t topArrayIndex = 0;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        YomoSequence *sequence = posCacheIt->second;
        sequence->setTopSegmentArrayIndex(topArrayIndex);
        inputIt = inputMap.find(sequence->getName());
        if (inputIt != inputMap.end()) {
            const Sequence::UpdateInfo *updateInfo = inputIt->second;
            newDimensions.push_back(*updateInfo);
        } else {
            currentIt = currentTopD.find(sequence->getName());
            assert(currentIt != currentTopD.end());
            newInfo._name = posCacheIt->first;
            newInfo._numSegments = currentIt->second;
            newDimensions.push_back(newInfo);
        }
        sequence->setNumTopSegments(newDimensions.back()._numSegments);
        topArrayIndex += newDimensions.back()._numSegments;
    }
    setGenomeTopDimensions(newDimensions);
}

void YomoGenome::updateBottomDimensions(const vector<Sequence::UpdateInfo> &bottomDimensions) {
    loadSequencePosCache();
    loadSequenceNameCache();
    vector<Sequence::UpdateInfo>::const_iterator i;
    map<string, YomoSequence *>::iterator cacheIt;
    map<string, const Sequence::UpdateInfo *> inputMap;
    map<string, hal_size_t> currentBottomD;
    // copy input into map, checking everything is already present
    for (i = bottomDimensions.begin(); i != bottomDimensions.end(); ++i) {
        const string &name = i->_name;
        cacheIt = _sequenceNameCache.find(name);
        if (cacheIt == _sequenceNameCache.end()) {
            throw hal_exception(string("Cannot update sequence ") + name + " because it is not present in "
                                                                           " genome " +
                                getName());
        }
        inputMap.insert(pair<string, const Sequence::UpdateInfo *>(name, &*i));
    }
    // keep a record of the number of segments in each existing
    // segment (these can get muddled as we add the new ones in the next
    // loop to be sure by getting them in one shot)
    // Note to self: iterating the map in this way skips zero-length
    // sequences (which are in the separate vector).  This is fine
    // here since we will never update them, but seems like it could be
    // dangerous if something were to change
    map<hal_size_t, YomoSequence *>::iterator posCacheIt;
    map<string, const Sequence::UpdateInfo *>::iterator inputIt;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        YomoSequence *sequence = posCacheIt->second;
        inputIt = inputMap.find(sequence->getName());
        if (inputIt == inputMap.end()) {
            currentBottomD.insert(pair<string, hal_size_t>(sequence->getName(), sequence->getNumBottomSegments()));
        }
    }
    // scan through existing sequences, updating as necessary
    // build summary of all new and unchanged dimensions in newDimensions
    // Note to self: iterating the map in this way skips zero-length
    // sequences (which are in the separate vector).  This is fine
    // here since we will never update them, but seems like it could be
    // dangerous if something were to change
    map<string, hal_size_t>::iterator currentIt;
    vector<Sequence::UpdateInfo> newDimensions;
    Sequence::UpdateInfo newInfo;
    hal_size_t bottomArrayIndex = 0;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        YomoSequence *sequence = posCacheIt->second;
        sequence->setBottomSegmentArrayIndex(bottomArrayIndex);
        inputIt = inputMap.find(sequence->getName());
        if (inputIt != inputMap.end()) {
            const Sequence::UpdateInfo *updateInfo = inputIt->second;
            newDimensions.push_back(*updateInfo);
        } else {
            currentIt = currentBottomD.find(sequence->getName());
            assert(currentIt != currentBottomD.end());
            newInfo._name = posCacheIt->first;
            newInfo._numSegments = currentIt->second;
            newDimensions.push_back(newInfo);
        }
        sequence->setNumBottomSegments(newDimensions.back()._numSegments);
        bottomArrayIndex += newDimensions.back()._numSegments;
    }
    setGenomeBottomDimensions(newDimensions);
}

void YomoGenome::setGenomeTopDimensions(const vector<Sequence::UpdateInfo> &topDimensions) {
    hal_size_t numTopSegments = 0;
    for (vector<Sequence::UpdateInfo>::const_iterator i = topDimensions.begin(); i != topDimensions.end(); ++i) {
        numTopSegments += i->_numSegments;
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(topArrayName);
        _group.unlink(topArrayName);
    } catch (H5::Exception &) {
    }
    _topArray.create(&_group, topArrayName, YomoTopSegment::dataType(), numTopSegments + 1, &_dcprops, _numChunksInArrayBuffer);
    reload();
}

void YomoGenome::setGenomeBottomDimensions(const vector<Sequence::UpdateInfo> &bottomDimensions) {
    hal_size_t numBottomSegments = 0;
    for (vector<Sequence::UpdateInfo>::const_iterator i = bottomDimensions.begin(); i != bottomDimensions.end(); ++i) {
        numBottomSegments += i->_numSegments;
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(bottomArrayName);
        _group.unlink(bottomArrayName);
    } catch (H5::Exception &) {
    }
    hal_size_t numChildren = _alignment->getChildNames(_name).size();

    // scale down the chunk size in order to keep chunks proportional to
    // the size of a bottom segment with two children.
    hsize_t chunk;
    _dcprops.getChunk(1, &chunk);
    double scale = numChildren < 10 ? 1. : 10. / numChildren;
    chunk *= scale;
    DSetCreatPropList botDC;
    botDC.copy(_dcprops);
    botDC.setChunk(1, &chunk);

    _bottomArray.create(&_group, bottomArrayName, YomoBottomSegment::dataType(numChildren), numBottomSegments + 1, &botDC,
                        _numChunksInArrayBuffer);
    reload();
}

hal_size_t YomoGenome::getNumSequences() const {
    assert(_sequenceIdxArray.getSize() == _sequenceNameArray.getSize() + 1);
    return _sequenceNameArray.getSize();
}

Sequence *YomoGenome::getSequence(const string &name) {
    loadSequenceNameCache();
    Sequence *sequence = NULL;
    map<string, YomoSequence *>::iterator mapIt = _sequenceNameCache.find(name);
    if (mapIt != _sequenceNameCache.end()) {
        sequence = mapIt->second;
    }
    return sequence;
}

const Sequence *YomoGenome::getSequence(const string &name) const {
    return const_cast<YomoGenome *>(this)->getSequence(name);
}

Sequence *YomoGenome::getSequenceBySite(hal_size_t position) {
    loadSequencePosCache();
    map<hal_size_t, YomoSequence *>::iterator i;
    i = _sequencePosCache.upper_bound(position);
    if (i != _sequencePosCache.end()) {
        if (position >= (hal_size_t)i->second->getStartPosition()) {
            assert(position < i->second->getStartPosition() + i->second->getSequenceLength());
            return i->second;
        }
    }
    return NULL;
}

const Sequence *YomoGenome::getSequenceBySite(hal_size_t position) const {
    return const_cast<YomoGenome *>(this)->getSequenceBySite(position);
}

SequenceIteratorPtr YomoGenome::getSequenceIterator(hal_index_t position) {
    assert(position <= (hal_index_t)_sequenceNameArray.getSize());
    YomoSequenceIterator *seqIt = new YomoSequenceIterator(this, position);
    return SequenceIteratorPtr(seqIt);
}

SequenceIteratorPtr YomoGenome::getSequenceIterator(hal_index_t position) const {
    return const_cast<YomoGenome *>(this)->getSequenceIterator(position);
}

MetaData *YomoGenome::getMetaData() {
    return _metaData;
}

const MetaData *YomoGenome::getMetaData() const {
    return _metaData;
}

bool YomoGenome::containsDNAArray() const {
    return _dnaArray.getSize() > 0;
}

const Alignment *YomoGenome::getAlignment() const {
    return _alignment;
}

Alignment *YomoGenome::getAlignment() {
    return _alignment;
}

// SEGMENTED SEQUENCE INTERFACE

const string &YomoGenome::getName() const {
    return _name;
}

hal_size_t YomoGenome::getSequenceLength() const {
    return _totalSequenceLength;
}

hal_size_t YomoGenome::getNumTopSegments() const {
    hal_size_t arraySize = _topArray.getSize();
    return arraySize > 0 ? arraySize - 1 : 0;
}

hal_size_t YomoGenome::getNumBottomSegments() const {
    hal_size_t arraySize = _bottomArray.getSize();
    return arraySize > 0 ? arraySize - 1 : 0;
}

TopSegmentIteratorPtr YomoGenome::getTopSegmentIterator(hal_index_t position) {
    assert(position <= (hal_index_t)getNumTopSegments());
    YomoTopSegment *topSeg = new YomoTopSegment(this, &_topArray, position);
    // ownership of topSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    TopSegmentIterator *topSegIt = new TopSegmentIterator(topSeg);
    return TopSegmentIteratorPtr(topSegIt);
}

TopSegmentIteratorPtr YomoGenome::getTopSegmentIterator(hal_index_t position) const {
    assert(position <= (hal_index_t)getNumTopSegments());
    YomoGenome *genome = const_cast<YomoGenome *>(this);
    YomoTopSegment *topSeg = new YomoTopSegment(genome, &genome->_topArray, position);
    // ownership of topSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    TopSegmentIterator *topSegIt = new TopSegmentIterator(topSeg);
    return TopSegmentIteratorPtr(topSegIt);
}

BottomSegmentIteratorPtr YomoGenome::getBottomSegmentIterator(hal_index_t position) {
    assert(position <= (hal_index_t)getNumBottomSegments());
    YomoBottomSegment *botSeg = new YomoBottomSegment(this, &_bottomArray, position);
    // ownership of botSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    BottomSegmentIterator *botSegIt = new BottomSegmentIterator(botSeg);
    return BottomSegmentIteratorPtr(botSegIt);
}

BottomSegmentIteratorPtr YomoGenome::getBottomSegmentIterator(hal_index_t position) const {
    assert(position <= (hal_index_t)getNumBottomSegments());
    YomoGenome *genome = const_cast<YomoGenome *>(this);
    YomoBottomSegment *botSeg = new YomoBottomSegment(genome, &genome->_bottomArray, position);
    // ownership of botSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    BottomSegmentIterator *botSegIt = new BottomSegmentIterator(botSeg);
    return BottomSegmentIteratorPtr(botSegIt);
}

DnaIteratorPtr YomoGenome::getDnaIterator(hal_index_t position) {
    assert(position / 2 <= (hal_index_t)_dnaArray.getSize());
    DnaIterator *dnaIt = new DnaIterator(this, _dnaAccess, position);
    return DnaIteratorPtr(dnaIt);
}

DnaIteratorPtr YomoGenome::getDnaIterator(hal_index_t position) const {
    return const_cast<YomoGenome *>(this)->getDnaIterator(position);
}

ColumnIteratorPtr YomoGenome::getColumnIterator(const set<const Genome *> *targets, hal_size_t maxInsertLength,
                                                hal_index_t position, hal_index_t lastPosition, bool noDupes, bool noAncestors,
                                                bool reverseStrand, bool unique, bool onlyOrthologs) const {
    hal_index_t lastIdx = lastPosition;
    if (lastPosition == NULL_INDEX) {
        lastIdx = (hal_index_t)(getSequenceLength() - 1);
    }
    if (position < 0 || lastPosition >= (hal_index_t)(getSequenceLength())) {
        throw hal_exception("YomoGenome::getColumnIterator: input indices (" + std::to_string(position) + ", " +
                            std::to_string(lastPosition) + ") out of bounds");
    }
    ColumnIterator *colIt = new ColumnIterator(this, targets, position, lastIdx, maxInsertLength, noDupes, noAncestors,
                                               reverseStrand, unique, onlyOrthologs);
    return ColumnIteratorPtr(colIt);
}

void YomoGenome::getString(string &outString) const {
    getSubString(outString, 0, getSequenceLength());
}

void YomoGenome::setString(const string &inString) {
    setSubString(inString, 0, getSequenceLength());
}

void YomoGenome::getSubString(string &outString, hal_size_t start, hal_size_t length) const {
    outString.resize(length);
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->readString(outString, length);
}

void YomoGenome::setSubString(const string &inString, hal_size_t start, hal_size_t length) {
    if (length != inString.length()) {
        throw hal_exception(string("setString: input string has different") + "length from target string in genome");
    }
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->writeString(inString, length);
}

RearrangementPtr YomoGenome::getRearrangement(hal_index_t position, hal_size_t gapLengthThreshold, double nThreshold,
                                              bool atomic) const {
    assert(position >= 0 && position < (hal_index_t)getNumTopSegments());
    TopSegmentIteratorPtr topIt = getTopSegmentIterator(position);
    Rearrangement *rea = new Rearrangement(this, gapLengthThreshold, nThreshold, atomic);
    rea->identifyFromLeftBreakpoint(topIt);
    return RearrangementPtr(rea);
}

GappedTopSegmentIteratorPtr YomoGenome::getGappedTopSegmentIterator(hal_index_t i, hal_size_t gapThreshold, bool atomic) const {
    TopSegmentIteratorPtr topSegIt = getTopSegmentIterator(i);
    GappedTopSegmentIterator *gapTopSegIt = new GappedTopSegmentIterator(topSegIt, gapThreshold, atomic);
    return GappedTopSegmentIteratorPtr(gapTopSegIt);
}

GappedBottomSegmentIteratorPtr YomoGenome::getGappedBottomSegmentIterator(hal_index_t i, hal_size_t childIdx,
                                                                          hal_size_t gapThreshold, bool atomic) const {
    BottomSegmentIteratorPtr botSegIt = getBottomSegmentIterator(i);
    GappedBottomSegmentIterator *gapBotSegIt = new GappedBottomSegmentIterator(botSegIt, childIdx, gapThreshold, atomic);
    return GappedBottomSegmentIteratorPtr(gapBotSegIt);
}

// LOCAL NON-INTERFACE METHODS

void YomoGenome::write() {
    _dnaArray.write();
    _topArray.write();
    _bottomArray.write();
    _metaData->write();
    _rup->write();
    _sequenceIdxArray.write();
    _sequenceNameArray.write();
}

void YomoGenome::read() {
    bool dnaLoaded = false;
    try {
        YOMODisableExceptionPrinting prDisable;
        _group.openDataSet(dnaArrayName);
        _dnaArray.load(&_group, dnaArrayName, _numChunksInArrayBuffer);
        dnaLoaded = true;
    } catch (H5::Exception &) {
    }

    try {
        YOMODisableExceptionPrinting prDisable;
        _group.openDataSet(topArrayName);
        _topArray.load(&_group, topArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        _group.openDataSet(bottomArrayName);
        _bottomArray.load(&_group, bottomArrayName, _numChunksInArrayBuffer);
        _numChildrenInBottomArray = YomoBottomSegment::numChildrenFromDataType(_bottomArray.getDataType());
    } catch (H5::Exception &) {
    }

    deleteSequenceCache();
    try {
        YOMODisableExceptionPrinting prDisable;
        _group.openDataSet(sequenceIdxArrayName);
        _sequenceIdxArray.load(&_group, sequenceIdxArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }
    try {
        YOMODisableExceptionPrinting prDisable;
        _group.openDataSet(sequenceNameArrayName);
        _sequenceNameArray.load(&_group, sequenceNameArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }

    readSequences();
    if (dnaLoaded) {
        _dnaAccess = DnaAccessPtr(new YOMODnaAccess(this, &_dnaArray, 0));
    }
}

void YomoGenome::readSequences() {
    deleteSequenceCache();
}

void YomoGenome::deleteSequenceCache() {
    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        map<hal_size_t, YomoSequence *>::iterator i;
        for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i) {
            delete i->second;
        }
        vector<YomoSequence *>::iterator z;
        for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z) {
            delete *z;
        }
    } else if (_sequenceNameCache.size() > 0) {
        map<string, YomoSequence *>::iterator i;
        for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i) {
            delete i->second;
        }
    }
    _sequencePosCache.clear();
    _zeroLenPosCache.clear();
    _sequenceNameCache.clear(); // I share my pointers with above.
}

void YomoGenome::loadSequencePosCache() const {
    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        return;
    }
    hal_size_t totalReadLen = 0;
    hal_size_t numSequences = _sequenceNameArray.getSize();

    if (_sequenceNameCache.size() > 0) {
        assert(_sequenceNameCache.size() == numSequences);
        map<std::string, YomoSequence *>::const_iterator i;
        for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i) {
            if (i->second->getSequenceLength() > 0) {
                _sequencePosCache.insert(pair<hal_size_t, YomoSequence *>(
                    i->second->getStartPosition() + i->second->getSequenceLength(), i->second));
                totalReadLen += i->second->getSequenceLength();
            } else {
                // FIXME: _zeroLenPosCache.push_back(i->second);
            }
        }
    } else {
        for (hal_size_t i = 0; i < numSequences; ++i) {
            YomoSequence *seq =
                new YomoSequence(const_cast<YomoGenome *>(this), const_cast<YomoExternalArray *>(&_sequenceIdxArray),
                                 const_cast<YomoExternalArray *>(&_sequenceNameArray), i);
            if (seq->getSequenceLength() > 0) {
                _sequencePosCache.insert(
                    pair<hal_size_t, YomoSequence *>(seq->getStartPosition() + seq->getSequenceLength(), seq));
                totalReadLen += seq->getSequenceLength();
            } else {
                // FIXME: _zeroLenPosCache.push_back(seq);
            }
        }
    }
    if (_totalSequenceLength > 0 && totalReadLen != _totalSequenceLength) {
        throw hal_exception("Sequences for genome " + getName() + " have total length " + std::to_string(totalReadLen) +
                            " but the (non-zero) DNA array contains " + std::to_string(_totalSequenceLength) +
                            " elements. This is an internal error " + "or the file is corrupt.");
    }
}

void YomoGenome::loadSequenceNameCache() const {
    if (_sequenceNameCache.size() > 0) {
        return;
    }
    hal_size_t numSequences = _sequenceNameArray.getSize();

    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        assert(_sequencePosCache.size() + _zeroLenPosCache.size() == numSequences);
        map<hal_size_t, YomoSequence *>::iterator i;
        for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i) {
            _sequenceNameCache.insert(pair<string, YomoSequence *>(i->second->getName(), i->second));
        }
        vector<YomoSequence *>::iterator z;
        for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z) {
            _sequenceNameCache.insert(pair<string, YomoSequence *>((*z)->getName(), (*z)));
        }
    } else {
        for (hal_size_t i = 0; i < numSequences; ++i) {
            YomoSequence *seq =
                new YomoSequence(const_cast<YomoGenome *>(this), const_cast<YomoExternalArray *>(&_sequenceIdxArray),
                                 const_cast<YomoExternalArray *>(&_sequenceNameArray), i);

            _sequenceNameCache.insert(pair<string, YomoSequence *>(seq->getName(), seq));
        }
    }
}

void YomoGenome::writeSequences(const vector<Sequence::Info> &sequenceDimensions) {
    deleteSequenceCache();
    vector<Sequence::Info>::const_iterator i;
    hal_size_t startPosition = 0;
    hal_size_t topArrayIndex = 0;
    hal_size_t bottomArrayIndex = 0;
    for (i = sequenceDimensions.begin(); i != sequenceDimensions.end(); ++i) {
        // Copy segment into YOMO array
        YomoSequence *seq = new YomoSequence(this, &_sequenceIdxArray, &_sequenceNameArray, i - sequenceDimensions.begin());
        // write all the Sequence::Info into the yomo sequence record
        seq->set(startPosition, *i, topArrayIndex, bottomArrayIndex);
        // Keep the object pointer in our caches
        if (seq->getSequenceLength() > 0) {
            _sequencePosCache.insert(pair<hal_size_t, YomoSequence *>(startPosition + i->_length, seq));
        } else {
            // FIXME: _zeroLenPosCache.push_back(seq);
        }
        _sequenceNameCache.insert(pair<string, YomoSequence *>(i->_name, seq));
        startPosition += i->_length;
        topArrayIndex += i->_numTopSegments;
        bottomArrayIndex += i->_numBottomSegments;
    }
}

void YomoGenome::resetBranchCaches() {
    _parentCache = NULL;
    _childCache.clear();
}

void YomoGenome::rename(const string &newName) {
    _group.move("/" + _name, "/" + newName);
    string newickStr = _alignment->getNewickTree();
    stTree *tree = stTree_parseNewickString(newickStr.c_str());
    stTree *node = stTree_findChild(tree, _name.c_str());
    stTree_setLabel(node, newName.c_str());
    _alignment->replaceNewickTree(stTree_getNewickTreeString(tree));
    stTree_destruct(tree);
}

void YomoGenome::renameSequence(const string &oldName, size_t index, const string &newName) {
    if (oldName.size() < newName.size()) {
        resizeNameArray(newName.size());
    }
    char *arrayBuffer = _sequenceNameArray.getUpdate(index);
    strcpy(arrayBuffer, newName.c_str());
    _sequenceNameArray.write();
    readSequences();
}

void YomoGenome::resizeNameArray(size_t newMaxSize) {
    size_t currentMaxSize = _sequenceNameArray.getDataType().getSize();
    if (newMaxSize > currentMaxSize) {
        vector<string> names;
        size_t numSequences = getNumSequences();
        for (hal_size_t i = 0; i < numSequences; ++i) {
            YomoSequence seq(const_cast<YomoGenome *>(this), const_cast<YomoExternalArray *>(&_sequenceIdxArray),
                             const_cast<YomoExternalArray *>(&_sequenceNameArray), i);
            names.push_back(seq.getName());
        }
        try {
            YOMODisableExceptionPrinting prDisable;
            DataSet d = _group.openDataSet(sequenceNameArrayName);
            _group.unlink(sequenceNameArrayName);
        } catch (H5::Exception &) {
        }

        _sequenceNameArray.create(&_group, sequenceNameArrayName, YomoSequence::nameDataType(newMaxSize + 1), numSequences,
                                  &_dcprops, _numChunksInArrayBuffer);
        for (size_t i = 0; i < numSequences; i++) {
            char *arrayBuffer = _sequenceNameArray.getUpdate(i);
            strcpy(arrayBuffer, names[i].c_str());
        }
    }
}
