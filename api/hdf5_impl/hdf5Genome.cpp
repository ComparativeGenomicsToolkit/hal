/* Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "hdf5Genome.h"
#include "H5Cpp.h"
#include "halBottomSegmentIterator.h"
#include "halColumnIterator.h"
#include "halDnaIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halRearrangement.h"
#include "halTopSegmentIterator.h"
#include "hdf5BottomSegment.h"
#include "hdf5Common.h"
#include "hdf5DnaDriver.h"
#include "hdf5Sequence.h"
#include "hdf5SequenceIterator.h"
#include "hdf5TopSegment.h"
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace hal;
using namespace std;
using namespace H5;

const string Hdf5Genome::dnaArrayName = "DNA_ARRAY";
const string Hdf5Genome::topArrayName = "TOP_ARRAY";
const string Hdf5Genome::bottomArrayName = "BOTTOM_ARRAY";
const string Hdf5Genome::sequenceIdxArrayName = "SEQIDX_ARRAY";
const string Hdf5Genome::sequenceNameArrayName = "SEQNAME_ARRAY";
const string Hdf5Genome::metaGroupName = "Meta";
const string Hdf5Genome::rupGroupName = "Rup";
const double Hdf5Genome::dnaChunkScale = 10.;
const hal_size_t Hdf5Genome::maxPosCache = 1000;
    
Hdf5Genome::Hdf5Genome(const string &name, Hdf5Alignment *alignment, PortableH5Location *h5Parent,
                       const DSetCreatPropList &dcProps, bool inMemory)
    : Genome(alignment, name), _alignment(alignment), _h5Parent(h5Parent), _name(name), _numChildrenInBottomArray(0),
      _totalSequenceLength(0), _numChunksInArrayBuffer(inMemory ? 0 : 1) {
    _dcprops.copy(dcProps);
    assert(!name.empty());
    assert(alignment != NULL && h5Parent != NULL);

    try {
        HDF5DisableExceptionPrinting prDisable;
        _group = h5Parent->openGroup(name);
    } catch (Exception &e) {
        _group = h5Parent->createGroup(name);
    }
    read();
    _metaData = new HDF5MetaData(&_group, metaGroupName);
    _rup = new HDF5MetaData(&_group, rupGroupName);

    _totalSequenceLength = _dnaArray.getSize() * 2;
    if (_totalSequenceLength > 0 && _rup->get(rupGroupName) == "1") {
        _totalSequenceLength -= 1;
    } else if (_totalSequenceLength == 0 && _sequenceIdxArray.getSize() > 0) {
        Hdf5Sequence lastSeq(this, &_sequenceIdxArray, &_sequenceNameArray, _sequenceNameArray.getSize() - 1);
        _totalSequenceLength = lastSeq.getEndPosition() + 1;
    }
}

Hdf5Genome::~Hdf5Genome() {
    delete _metaData;
    delete _rup;
    deleteSequenceCache();
}

// GENOME INTERFACE

void Hdf5Genome::setDimensions(const vector<Sequence::Info> &sequenceDimensions, bool storeDNAArrays) {
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
        HDF5DisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(dnaArrayName);
        _group.unlink(dnaArrayName);
    } catch (H5::Exception &) {
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(sequenceIdxArrayName);
        _group.unlink(sequenceIdxArrayName);
    } catch (H5::Exception &) {
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
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
        _dnaAccess = DnaAccessPtr(new HDF5DnaAccess(this, &_dnaArray, 0));
    }
    if (totalSeq > 0) {
        _sequenceIdxArray.create(&_group, sequenceIdxArrayName, Hdf5Sequence::idxDataType(), totalSeq + 1, &_dcprops,
                                 _numChunksInArrayBuffer);

        _sequenceNameArray.create(&_group, sequenceNameArrayName, Hdf5Sequence::nameDataType(maxName + 1), totalSeq, &_dcprops,
                                  _numChunksInArrayBuffer);

        writeSequences(sequenceDimensions);
    }

    // Do the same as above for the segments.
    setGenomeTopDimensions(topDimensions);
    setGenomeBottomDimensions(bottomDimensions);

    reload();
}

void Hdf5Genome::updateTopDimensions(const vector<Sequence::UpdateInfo> &topDimensions) {
    loadSequencePosCache();
    loadSequenceNameCache();
    vector<Sequence::UpdateInfo>::const_iterator i;
    map<string, Hdf5Sequence *>::iterator cacheIt;
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
    map<hal_size_t, Hdf5Sequence *>::iterator posCacheIt;
    map<string, const Sequence::UpdateInfo *>::iterator inputIt;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        Hdf5Sequence *sequence = posCacheIt->second;
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
        Hdf5Sequence *sequence = posCacheIt->second;
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

void Hdf5Genome::updateBottomDimensions(const vector<Sequence::UpdateInfo> &bottomDimensions) {
    loadSequencePosCache();
    loadSequenceNameCache();
    vector<Sequence::UpdateInfo>::const_iterator i;
    map<string, Hdf5Sequence *>::iterator cacheIt;
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
    map<hal_size_t, Hdf5Sequence *>::iterator posCacheIt;
    map<string, const Sequence::UpdateInfo *>::iterator inputIt;
    for (posCacheIt = _sequencePosCache.begin(); posCacheIt != _sequencePosCache.end(); ++posCacheIt) {
        Hdf5Sequence *sequence = posCacheIt->second;
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
        Hdf5Sequence *sequence = posCacheIt->second;
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

void Hdf5Genome::setGenomeTopDimensions(const vector<Sequence::UpdateInfo> &topDimensions) {
    hal_size_t numTopSegments = 0;
    for (vector<Sequence::UpdateInfo>::const_iterator i = topDimensions.begin(); i != topDimensions.end(); ++i) {
        numTopSegments += i->_numSegments;
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
        DataSet d = _group.openDataSet(topArrayName);
        _group.unlink(topArrayName);
    } catch (H5::Exception &) {
    }
    _topArray.create(&_group, topArrayName, Hdf5TopSegment::dataType(), numTopSegments + 1, &_dcprops, _numChunksInArrayBuffer);
    reload();
}

void Hdf5Genome::setGenomeBottomDimensions(const vector<Sequence::UpdateInfo> &bottomDimensions) {
    hal_size_t numBottomSegments = 0;
    for (vector<Sequence::UpdateInfo>::const_iterator i = bottomDimensions.begin(); i != bottomDimensions.end(); ++i) {
        numBottomSegments += i->_numSegments;
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
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

    _bottomArray.create(&_group, bottomArrayName, Hdf5BottomSegment::dataType(numChildren), numBottomSegments + 1, &botDC,
                        _numChunksInArrayBuffer);
    reload();
}

hal_size_t Hdf5Genome::getNumSequences() const {
    assert(_sequenceIdxArray.getSize() == _sequenceNameArray.getSize() + 1);
    return _sequenceNameArray.getSize();
}

Sequence *Hdf5Genome::getSequence(const string &name) {
    loadSequenceNameCache();
    Sequence *sequence = NULL;
    map<string, Hdf5Sequence *>::iterator mapIt = _sequenceNameCache.find(name);
    if (mapIt != _sequenceNameCache.end()) {
        sequence = mapIt->second;
    }
    return sequence;
}

const Sequence *Hdf5Genome::getSequence(const string &name) const {
    return const_cast<Hdf5Genome *>(this)->getSequence(name);
}

Sequence *Hdf5Genome::getSequenceBySite(hal_size_t position) {
    hal_size_t numSequences = _sequenceNameArray.getSize();
    if (numSequences <= maxPosCache) {
        loadSequencePosCache();
    }
    
    // check the cache
    map<hal_size_t, Hdf5Sequence *>::iterator i;
    i = _sequencePosCache.upper_bound(position);
    if (i != _sequencePosCache.end()) {
        if (position >= (hal_size_t)i->second->getStartPosition() && position <= (hal_size_t)i->second->getEndPosition()) {
            return i->second;
        }
    }
    
    // it wasn't in the pos cache.  we dig it out of the array using binary search
    hal_index_t low_idx = 0;
    hal_index_t high_idx = (hal_index_t)numSequences - 1;
    
    while (low_idx <= high_idx) {
        hal_index_t mid_idx = (low_idx + high_idx) / 2;
        Hdf5Sequence *seq =
            new Hdf5Sequence(const_cast<Hdf5Genome *>(this), const_cast<Hdf5ExternalArray *>(&_sequenceIdxArray),
                             const_cast<Hdf5ExternalArray *>(&_sequenceNameArray), mid_idx);
        if (position >= seq->getStartPosition() && position <= seq->getEndPosition()) {
            _sequencePosCache[seq->getStartPosition() + seq->getSequenceLength()] = seq;
            return seq;
        } else {
            bool less = position < seq->getStartPosition();
            delete seq;
            if (less) {
                high_idx = mid_idx - 1;
            } else {
                low_idx = mid_idx + 1;
            }
        }
    }

    // something went wrong (sequences stored out of order?!), fall back on the cache
    if (numSequences > maxPosCache) {
        cerr << "Binary search failed on " << getName() << ":" << position << " reverting to full cache" << endl;
        loadSequencePosCache();
        i = _sequencePosCache.upper_bound(position);
        if (i != _sequencePosCache.end()) {
            if (position >= (hal_size_t)i->second->getStartPosition() && position <= (hal_size_t)i->second->getEndPosition()) {
                return i->second;
            }
        }
    }
    throw hal_exception("Position cache fail on " + getName() + ":" + std::to_string(position));
    return NULL;
}

const Sequence *Hdf5Genome::getSequenceBySite(hal_size_t position) const {
    return const_cast<Hdf5Genome *>(this)->getSequenceBySite(position);
}

SequenceIteratorPtr Hdf5Genome::getSequenceIterator(hal_index_t position) {
    assert(position <= (hal_index_t)_sequenceNameArray.getSize());
    Hdf5SequenceIterator *seqIt = new Hdf5SequenceIterator(this, position);
    return SequenceIteratorPtr(seqIt);
}

SequenceIteratorPtr Hdf5Genome::getSequenceIterator(hal_index_t position) const {
    return const_cast<Hdf5Genome *>(this)->getSequenceIterator(position);
}

MetaData *Hdf5Genome::getMetaData() {
    return _metaData;
}

const MetaData *Hdf5Genome::getMetaData() const {
    return _metaData;
}

bool Hdf5Genome::containsDNAArray() const {
    return _dnaArray.getSize() > 0;
}

const Alignment *Hdf5Genome::getAlignment() const {
    return _alignment;
}

Alignment *Hdf5Genome::getAlignment() {
    return _alignment;
}

// SEGMENTED SEQUENCE INTERFACE

const string &Hdf5Genome::getName() const {
    return _name;
}

hal_size_t Hdf5Genome::getSequenceLength() const {
    return _totalSequenceLength;
}

hal_size_t Hdf5Genome::getNumTopSegments() const {
    hal_size_t arraySize = _topArray.getSize();
    return arraySize > 0 ? arraySize - 1 : 0;
}

hal_size_t Hdf5Genome::getNumBottomSegments() const {
    hal_size_t arraySize = _bottomArray.getSize();
    return arraySize > 0 ? arraySize - 1 : 0;
}

TopSegmentIteratorPtr Hdf5Genome::getTopSegmentIterator(hal_index_t position) {
    assert(position <= (hal_index_t)getNumTopSegments());
    Hdf5TopSegment *topSeg = new Hdf5TopSegment(this, &_topArray, position);
    // ownership of topSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    TopSegmentIterator *topSegIt = new TopSegmentIterator(topSeg);
    return TopSegmentIteratorPtr(topSegIt);
}

TopSegmentIteratorPtr Hdf5Genome::getTopSegmentIterator(hal_index_t position) const {
    assert(position <= (hal_index_t)getNumTopSegments());
    Hdf5Genome *genome = const_cast<Hdf5Genome *>(this);
    Hdf5TopSegment *topSeg = new Hdf5TopSegment(genome, &genome->_topArray, position);
    // ownership of topSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    TopSegmentIterator *topSegIt = new TopSegmentIterator(topSeg);
    return TopSegmentIteratorPtr(topSegIt);
}

BottomSegmentIteratorPtr Hdf5Genome::getBottomSegmentIterator(hal_index_t position) {
    assert(position <= (hal_index_t)getNumBottomSegments());
    Hdf5BottomSegment *botSeg = new Hdf5BottomSegment(this, &_bottomArray, position);
    // ownership of botSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    BottomSegmentIterator *botSegIt = new BottomSegmentIterator(botSeg);
    return BottomSegmentIteratorPtr(botSegIt);
}

BottomSegmentIteratorPtr Hdf5Genome::getBottomSegmentIterator(hal_index_t position) const {
    assert(position <= (hal_index_t)getNumBottomSegments());
    Hdf5Genome *genome = const_cast<Hdf5Genome *>(this);
    Hdf5BottomSegment *botSeg = new Hdf5BottomSegment(genome, &genome->_bottomArray, position);
    // ownership of botSeg is passed into newIt, whose lifespan is
    // governed by the returned smart pointer
    BottomSegmentIterator *botSegIt = new BottomSegmentIterator(botSeg);
    return BottomSegmentIteratorPtr(botSegIt);
}

DnaIteratorPtr Hdf5Genome::getDnaIterator(hal_index_t position) {
    assert(position / 2 <= (hal_index_t)_dnaArray.getSize());
    DnaIterator *dnaIt = new DnaIterator(this, _dnaAccess, position);
    return DnaIteratorPtr(dnaIt);
}

DnaIteratorPtr Hdf5Genome::getDnaIterator(hal_index_t position) const {
    return const_cast<Hdf5Genome *>(this)->getDnaIterator(position);
}

ColumnIteratorPtr Hdf5Genome::getColumnIterator(const set<const Genome *> *targets, hal_size_t maxInsertLength,
                                                hal_index_t position, hal_index_t lastPosition, bool noDupes, bool noAncestors,
                                                bool reverseStrand, bool unique, bool onlyOrthologs) const {
    hal_index_t lastIdx = lastPosition;
    if (lastPosition == NULL_INDEX) {
        lastIdx = (hal_index_t)(getSequenceLength() - 1);
    }
    if (position < 0 || lastPosition >= (hal_index_t)(getSequenceLength())) {
        throw hal_exception("Hdf5Genome::getColumnIterator: input indices (" + std::to_string(position) + ", " +
                            std::to_string(lastPosition) + ") out of bounds");
    }
    ColumnIterator *colIt = new ColumnIterator(this, targets, position, lastIdx, maxInsertLength, noDupes, noAncestors,
                                               reverseStrand, unique, onlyOrthologs);
    return ColumnIteratorPtr(colIt);
}

void Hdf5Genome::getString(string &outString) const {
    getSubString(outString, 0, getSequenceLength());
}

void Hdf5Genome::setString(const string &inString) {
    setSubString(inString, 0, getSequenceLength());
}

void Hdf5Genome::getSubString(string &outString, hal_size_t start, hal_size_t length) const {
    outString.resize(length);
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->readString(outString, length);
}

void Hdf5Genome::setSubString(const string &inString, hal_size_t start, hal_size_t length) {
    if (length != inString.length()) {
        throw hal_exception(string("setString: input string has different") + "length from target string in genome");
    }
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->writeString(inString, length);
}

RearrangementPtr Hdf5Genome::getRearrangement(hal_index_t position, hal_size_t gapLengthThreshold, double nThreshold,
                                              bool atomic) const {
    assert(position >= 0 && position < (hal_index_t)getNumTopSegments());
    TopSegmentIteratorPtr topIt = getTopSegmentIterator(position);
    Rearrangement *rea = new Rearrangement(this, gapLengthThreshold, nThreshold, atomic);
    rea->identifyFromLeftBreakpoint(topIt);
    return RearrangementPtr(rea);
}

GappedTopSegmentIteratorPtr Hdf5Genome::getGappedTopSegmentIterator(hal_index_t i, hal_size_t gapThreshold, bool atomic) const {
    TopSegmentIteratorPtr topSegIt = getTopSegmentIterator(i);
    GappedTopSegmentIterator *gapTopSegIt = new GappedTopSegmentIterator(topSegIt, gapThreshold, atomic);
    return GappedTopSegmentIteratorPtr(gapTopSegIt);
}

GappedBottomSegmentIteratorPtr Hdf5Genome::getGappedBottomSegmentIterator(hal_index_t i, hal_size_t childIdx,
                                                                          hal_size_t gapThreshold, bool atomic) const {
    BottomSegmentIteratorPtr botSegIt = getBottomSegmentIterator(i);
    GappedBottomSegmentIterator *gapBotSegIt = new GappedBottomSegmentIterator(botSegIt, childIdx, gapThreshold, atomic);
    return GappedBottomSegmentIteratorPtr(gapBotSegIt);
}

// LOCAL NON-INTERFACE METHODS

void Hdf5Genome::write() {
    _dnaArray.write();
    _topArray.write();
    _bottomArray.write();
    _metaData->write();
    _rup->write();
    _sequenceIdxArray.write();
    _sequenceNameArray.write();
}

void Hdf5Genome::read() {
    bool dnaLoaded = false;
    try {
        HDF5DisableExceptionPrinting prDisable;
        _group.openDataSet(dnaArrayName);
        _dnaArray.load(&_group, dnaArrayName, _numChunksInArrayBuffer);
        dnaLoaded = true;
    } catch (H5::Exception &) {
    }

    try {
        HDF5DisableExceptionPrinting prDisable;
        _group.openDataSet(topArrayName);
        _topArray.load(&_group, topArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
        _group.openDataSet(bottomArrayName);
        _bottomArray.load(&_group, bottomArrayName, _numChunksInArrayBuffer);
        _numChildrenInBottomArray = Hdf5BottomSegment::numChildrenFromDataType(_bottomArray.getDataType());
    } catch (H5::Exception &) {
    }

    deleteSequenceCache();
    try {
        HDF5DisableExceptionPrinting prDisable;
        _group.openDataSet(sequenceIdxArrayName);
        _sequenceIdxArray.load(&_group, sequenceIdxArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }
    try {
        HDF5DisableExceptionPrinting prDisable;
        _group.openDataSet(sequenceNameArrayName);
        _sequenceNameArray.load(&_group, sequenceNameArrayName, _numChunksInArrayBuffer);
    } catch (H5::Exception &) {
    }

    readSequences();
    if (dnaLoaded) {
        _dnaAccess = DnaAccessPtr(new HDF5DnaAccess(this, &_dnaArray, 0));
    }
}

void Hdf5Genome::readSequences() {
    deleteSequenceCache();
}

void Hdf5Genome::deleteSequenceCache() {
    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        map<hal_size_t, Hdf5Sequence *>::iterator i;
        for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i) {
            delete i->second;
        }
        vector<Hdf5Sequence *>::iterator z;
        for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z) {
            delete *z;
        }
    } else if (_sequenceNameCache.size() > 0) {
        map<string, Hdf5Sequence *>::iterator i;
        for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i) {
            delete i->second;
        }
    }
    _sequencePosCache.clear();
    _zeroLenPosCache.clear();
    _sequenceNameCache.clear(); // I share my pointers with above.
}

void Hdf5Genome::loadSequencePosCache() const {
    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        return;
    }
    hal_size_t totalReadLen = 0;
    hal_size_t numSequences = _sequenceNameArray.getSize();

    if (_sequenceNameCache.size() > 0) {
        assert(_sequenceNameCache.size() == numSequences);
        map<std::string, Hdf5Sequence *>::const_iterator i;
        for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i) {
            if (i->second->getSequenceLength() > 0) {
                _sequencePosCache.insert(pair<hal_size_t, Hdf5Sequence *>(
                    i->second->getStartPosition() + i->second->getSequenceLength(), i->second));
                totalReadLen += i->second->getSequenceLength();
            } else {
                // FIXME: _zeroLenPosCache.push_back(i->second);
            }
        }
    } else {
        for (hal_size_t i = 0; i < numSequences; ++i) {
            Hdf5Sequence *seq =
                new Hdf5Sequence(const_cast<Hdf5Genome *>(this), const_cast<Hdf5ExternalArray *>(&_sequenceIdxArray),
                                 const_cast<Hdf5ExternalArray *>(&_sequenceNameArray), i);
            if (seq->getSequenceLength() > 0) {
                _sequencePosCache.insert(
                    pair<hal_size_t, Hdf5Sequence *>(seq->getStartPosition() + seq->getSequenceLength(), seq));
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

void Hdf5Genome::loadSequenceNameCache() const {
    if (_sequenceNameCache.size() > 0) {
        return;
    }
    hal_size_t numSequences = _sequenceNameArray.getSize();

    if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0) {
        assert(_sequencePosCache.size() + _zeroLenPosCache.size() == numSequences);
        map<hal_size_t, Hdf5Sequence *>::iterator i;
        for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i) {
            _sequenceNameCache.insert(pair<string, Hdf5Sequence *>(i->second->getName(), i->second));
        }
        vector<Hdf5Sequence *>::iterator z;
        for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z) {
            _sequenceNameCache.insert(pair<string, Hdf5Sequence *>((*z)->getName(), (*z)));
        }
    } else {
        for (hal_size_t i = 0; i < numSequences; ++i) {
            Hdf5Sequence *seq =
                new Hdf5Sequence(const_cast<Hdf5Genome *>(this), const_cast<Hdf5ExternalArray *>(&_sequenceIdxArray),
                                 const_cast<Hdf5ExternalArray *>(&_sequenceNameArray), i);

            _sequenceNameCache.insert(pair<string, Hdf5Sequence *>(seq->getName(), seq));
        }
    }
}

void Hdf5Genome::writeSequences(const vector<Sequence::Info> &sequenceDimensions) {
    deleteSequenceCache();
    vector<Sequence::Info>::const_iterator i;
    hal_size_t startPosition = 0;
    hal_size_t topArrayIndex = 0;
    hal_size_t bottomArrayIndex = 0;
    for (i = sequenceDimensions.begin(); i != sequenceDimensions.end(); ++i) {
        // Copy segment into HDF5 array
        Hdf5Sequence *seq = new Hdf5Sequence(this, &_sequenceIdxArray, &_sequenceNameArray, i - sequenceDimensions.begin());
        // write all the Sequence::Info into the hdf5 sequence record
        seq->set(startPosition, *i, topArrayIndex, bottomArrayIndex);
        // Keep the object pointer in our caches
        if (seq->getSequenceLength() > 0) {
            _sequencePosCache.insert(pair<hal_size_t, Hdf5Sequence *>(startPosition + i->_length, seq));
        } else {
            // FIXME: _zeroLenPosCache.push_back(seq);
        }
        _sequenceNameCache.insert(pair<string, Hdf5Sequence *>(i->_name, seq));
        startPosition += i->_length;
        topArrayIndex += i->_numTopSegments;
        bottomArrayIndex += i->_numBottomSegments;
    }
}

void Hdf5Genome::resetBranchCaches() {
    _parentCache = NULL;
    _childCache.clear();
}

void Hdf5Genome::rename(const string &newName) {
    _group.move("/" + _name, "/" + newName);
    string newickStr = _alignment->getNewickTree();
    stTree *tree = stTree_parseNewickString(newickStr.c_str());
    stTree *node = stTree_findChild(tree, _name.c_str());
    stTree_setLabel(node, newName.c_str());
    _alignment->replaceNewickTree(stTree_getNewickTreeString(tree));
    stTree_destruct(tree);
}

void Hdf5Genome::renameSequence(const string &oldName, size_t index, const string &newName) {
    if (oldName.size() < newName.size()) {
        resizeNameArray(newName.size() + 1);
    }
    char *arrayBuffer = _sequenceNameArray.getUpdate(index);
    strcpy(arrayBuffer, newName.c_str());
    _sequenceNameArray.write();
    readSequences();
}

void Hdf5Genome::resizeNameArray(size_t newMaxSize) {
    size_t currentMaxSize = _sequenceNameArray.getDataType().getSize();
    if (newMaxSize > currentMaxSize) {
        vector<string> names;
        size_t numSequences = getNumSequences();
        for (hal_size_t i = 0; i < numSequences; ++i) {
            Hdf5Sequence seq(const_cast<Hdf5Genome *>(this), const_cast<Hdf5ExternalArray *>(&_sequenceIdxArray),
                             const_cast<Hdf5ExternalArray *>(&_sequenceNameArray), i);
            names.push_back(seq.getName());
        }
        try {
            HDF5DisableExceptionPrinting prDisable;
            DataSet d = _group.openDataSet(sequenceNameArrayName);
            _group.unlink(sequenceNameArrayName);
        } catch (H5::Exception &) {
        }

        _sequenceNameArray.create(&_group, sequenceNameArrayName, Hdf5Sequence::nameDataType(newMaxSize), numSequences,
                                  &_dcprops, _numChunksInArrayBuffer);
        for (size_t i = 0; i < numSequences; i++) {
            char *arrayBuffer = _sequenceNameArray.getUpdate(i);
            strcpy(arrayBuffer, names[i].c_str());
        }
    }
}
