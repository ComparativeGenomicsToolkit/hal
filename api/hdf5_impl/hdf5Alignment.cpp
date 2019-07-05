/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hdf5Alignment.h"
#include "halCLParser.h"
#include "halCommon.h"
#include "halSequenceIterator.h"
#include "hdf5Common.h"
#include "hdf5Genome.h"
#include "hdf5MetaData.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
extern "C" {
#include "sonLibTree.h"
}

using namespace hal;
using namespace std;
using namespace H5;

#ifdef ENABLE_UDC
#include "hdf5UDCFuseDriver.h"
static const hid_t UDC_FUSE_DRIVER_ID = H5FD_udc_fuse_init();
#endif

/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /Meta */
const H5std_string Hdf5Alignment::MetaGroupName = "Meta";
const H5std_string Hdf5Alignment::TreeGroupName = "Phylogeny";
const H5std_string Hdf5Alignment::GenomesGroupName = "Genomes";
const H5std_string Hdf5Alignment::VersionGroupName = "Verison";

const hsize_t Hdf5Alignment::DefaultChunkSize = 1000;
const hsize_t Hdf5Alignment::DefaultCompression = 2;
const hsize_t Hdf5Alignment::DefaultCacheMDCElems = 113;
const hsize_t Hdf5Alignment::DefaultCacheRDCElems = 599999;
const hsize_t Hdf5Alignment::DefaultCacheRDCBytes = 15728640;
const double Hdf5Alignment::DefaultCacheW0 = 0.75;
const bool Hdf5Alignment::DefaultInMemory = false;

/* check if first bit of file has HDF5 header */
bool hal::Hdf5Alignment::isHdf5File(const std::string &initialBytes) {
    static const std::string HDF5_MAGIC = "\x89HDF\r";
    return initialBytes.compare(0, HDF5_MAGIC.size(), HDF5_MAGIC) == 0;
}

/* construction default flags */
static int hdf5DefaultFlags(unsigned mode) {
    if (mode & CREATE_ACCESS) {
        return H5F_ACC_TRUNC;
    } else if (mode & WRITE_ACCESS) {
        return H5F_ACC_RDWR;
    } else {
        return H5F_ACC_RDONLY;
    }
}

Hdf5Alignment::Hdf5Alignment(const string &alignmentPath, unsigned mode, const H5::FileCreatPropList &fileCreateProps,
                             const H5::FileAccPropList &fileAccessProps, const H5::DSetCreatPropList &datasetCreateProps,
                             bool inMemory)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(inMemory), _metaData(NULL), _tree(NULL), _dirty(false) {
    _cprops.copy(fileCreateProps);
    _aprops.copy(fileAccessProps);
    _dcprops.copy(datasetCreateProps);
    if (_inMemory) {
        setInMemory();
    }
    if (_mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

Hdf5Alignment::Hdf5Alignment(const std::string &alignmentPath, unsigned mode, const CLParser *parser)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _file(NULL), _flags(hdf5DefaultFlags(_mode)),
      _inMemory(false), _metaData(NULL), _tree(NULL), _dirty(false) {
    initializeFromOptions(parser);
    if (_inMemory) {
        setInMemory();
    }
    if (_mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

Hdf5Alignment::~Hdf5Alignment() {
    close();
}

void Hdf5Alignment::defineOptions(CLParser *parser, unsigned mode) {
    if ((mode & CREATE_ACCESS) || (mode & WRITE_ACCESS)) {
        parser->addOption("hdf5Chunk", "hdf5 chunk size", DefaultChunkSize);
        parser->addOption("chunk", "obsolete name for --hdf5Chunk ", DefaultChunkSize);

        parser->addOption("hdf5Compression", "hdf5 compression factor [0:none - 9:max]", DefaultCompression);
        parser->addOption("deflate", "obsolete name for --hdf5Compression", DefaultCompression);
    }
    parser->addOption("hdf5CacheMDC", "number of metadata slots in hdf5 cache", DefaultCacheMDCElems);
    parser->addOption("cacheMDC", "obsolete name for --hdf5CacheMDC ", DefaultCacheMDCElems);

    parser->addOption("hdf5CacheRDC", "number of regular slots in hdf5 cache.  should be"
                                      " a prime number ~= 10 * DefaultCacheRDCBytes / chunk",
                      DefaultCacheRDCElems);
    parser->addOption("cacheRDC", "obsolete name for --hdf5CacheRDC", Hdf5Alignment::DefaultCacheRDCElems);

    parser->addOption("hdf5CacheBytes", "maximum size in bytes of regular hdf5 cache", DefaultCacheRDCBytes);
    parser->addOption("cacheBytes", "obsolete name for --hdf5CacheBytes", DefaultCacheRDCBytes);

    parser->addOption("hdf5CacheW0", "w0 parameter for hdf5 cache", DefaultCacheW0);
    parser->addOption("cacheW0", "obsolete name for --hdf5CacheW0", DefaultCacheW0);

    parser->addOptionFlag("hdf5InMemory", "load all data in memory (and disable hdf5 cache)", DefaultInMemory);
    parser->addOptionFlag("inMemory", "obsolete name for --hdf5InMemory", DefaultInMemory);
}

/* initialize class from options */
void Hdf5Alignment::initializeFromOptions(const CLParser *parser) {
    _cprops.copy(H5::FileCreatPropList::DEFAULT);
    _dcprops.copy(H5::DSetCreatPropList::DEFAULT);
    _aprops.copy(H5::FileAccPropList::DEFAULT);
    _inMemory = parser->getFlagAlt("hdf5InMemory", "inMemory");
    if ((_mode & CREATE_ACCESS) || (_mode & WRITE_ACCESS)) {
        // these are only available on create
        hsize_t chunk = parser->getOptionAlt<hsize_t>("hdf5Chunk", "chunk");
        _dcprops.setChunk(1, &chunk);
        _dcprops.setDeflate(parser->getOptionAlt<hsize_t>("hdf5Compression", "deflate"));
    }
    _aprops.setCache(
        parser->getOptionAlt<hsize_t>("hdf5CacheMDC", "cacheMDC"), parser->getOptionAlt<hsize_t>("hdf5CacheRDC", "cacheRDC"),
        parser->getOptionAlt<hsize_t>("hdf5CacheBytes", "cacheBytes"), parser->getOptionAlt<double>("hdf5CacheW0", "cacheW0"));
}

/* set properties for in-memory access */
void Hdf5Alignment::setInMemory() {
    int mdc;
    size_t rdc, rdcb;
    double w0;
    _aprops.getCache(mdc, rdc, rdcb, w0);
    _aprops.setCache(mdc, 0, 0, 0.0);
}

void Hdf5Alignment::create() {
    if (not ofstream(_alignmentPath.c_str())) { // FIXME report errno
        throw hal_exception("Unable to open " + _alignmentPath);
    }

    _file = new H5File(_alignmentPath.c_str(), _flags, _cprops, _aprops);
    _file->createGroup(MetaGroupName);
    _file->createGroup(TreeGroupName);
    _file->createGroup(GenomesGroupName);
    _file->createGroup(VersionGroupName);
    _metaData = new HDF5MetaData(_file, MetaGroupName);
    _tree = NULL;
    _dirty = true;
    writeVersion();
}

void Hdf5Alignment::open() {
#ifdef ENABLE_UDC
    // check for url rather than UDC cache, since default UDC cache directory could be used.
    if (((_mode & WRITE_ACCESS) == 0) and isUrl(_alignmentPath)) {
        _aprops.setDriver(UDC_FUSE_DRIVER_ID, NULL);
    }
#else
    if (not ifstream(_alignmentPath.c_str())) { // FIXME report errno
        throw hal_errno_exception("Unable to open " + _alignmentPath, errno);
    }
#endif
    _file = new H5File(_alignmentPath.c_str(), _flags, _cprops, _aprops);
    if (!compatibleWithVersion(getVersion())) {
        throw hal_exception("HAL API v" + HAL_VERSION + " incompatible with format v" + getVersion() + " HAL file.");
    }
    _metaData = new HDF5MetaData(_file, MetaGroupName);
    loadTree();
    if (stTree_getNumNodes(_tree) > 199) {
        // We disable the HDF5 chunk cache with large numbers of species
        // (over 200 genomes, or 100 leaves with a binary tree), because
        // every genome opened will cause several new chunk caches to be
        // initialized. This causes massive memory usage (well over 17
        // GB for just running halStats on a 250-genome alignment).
        close();
        delete _file;
        _aprops.setCache(0, 0, 0, 0.);
        _file = new H5File(_alignmentPath.c_str(), _flags, _cprops, _aprops);
        delete _metaData;
        _metaData = new HDF5MetaData(_file, MetaGroupName);
        loadTree();
    }
}

void Hdf5Alignment::close() {
    if (_file != NULL) {
        if (not isReadOnly()) {
            writeTree();
        }
        if (_tree != NULL) {
            stTree_destruct(_tree);
            _tree = NULL;
        }
        // todo: make sure there's no memory leak with metadata
        // smart pointer should prevent
        if (_metaData != NULL) {
            if (not isReadOnly()) {
                _metaData->write();
            }
            delete _metaData;
            _metaData = NULL;
        }
        writeVersion();
        map<string, Hdf5Genome *>::iterator mapIt;
        for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt) {
            Hdf5Genome *genome = mapIt->second;
            if (not isReadOnly()) {
                genome->write();
            }
            delete genome;
        }
        _openGenomes.clear();
        if (not isReadOnly()) {
            _file->flush(H5F_SCOPE_LOCAL);
        }
        _file->close();
        delete _file;
        _file = NULL;
    } else {
        assert(_tree == NULL);
        assert(_openGenomes.empty() == true);
    }
}

Genome *Hdf5Alignment::insertGenome(const string &name, const string &parentName, const string &childName,
                                    double upperBranchLength) {
    if (name.empty() == true || parentName.empty() || childName.empty()) {
        throw hal_exception("name can't be empty");
    }
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt != _nodeMap.end()) {
        throw hal_exception("node " + name + " already exists");
    }
    findIt = _nodeMap.find(childName);
    if (findIt == _nodeMap.end()) {
        throw hal_exception("child " + childName + " not found in tree");
    }
    stTree *child = findIt->second;
    stTree *parent = stTree_getParent(child);
    if (stTree_getLabel(parent) != parentName) {
        throw hal_exception("no edge between " + parentName + " and " + childName);
    }
    double existingBranchLength = getBranchLength(parentName, childName);
    stTree *newNode = stTree_construct();
    stTree_setLabel(newNode, name.c_str());
    stTree_setParent(newNode, parent);
    stTree_setBranchLength(newNode, upperBranchLength);
    _nodeMap.insert(pair<string, stTree *>(name, newNode));
    double lowerBranchLength = existingBranchLength - upperBranchLength;
    stTree_setParent(child, newNode);
    stTree_setBranchLength(child, lowerBranchLength);

    Hdf5Genome *genome = new Hdf5Genome(name, this, _file, _dcprops, _inMemory);
    _openGenomes.insert(pair<string, Hdf5Genome *>(name, genome));
    _dirty = true;
    return genome;
}

Genome *Hdf5Alignment::addLeafGenome(const string &name, const string &parentName, double branchLength) {
    if (name.empty() == true || parentName.empty()) {
        throw hal_exception("name can't be empty");
    }
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt != _nodeMap.end()) {
        throw hal_exception(string("node ") + name + " already exists");
    }
    findIt = _nodeMap.find(parentName);
    if (findIt == _nodeMap.end()) {
        throw hal_exception(string("parent ") + parentName + " not found in tree");
    }
    stTree *parent = findIt->second;
    stTree *childNode = stTree_construct();
    stTree_setLabel(childNode, name.c_str());
    stTree_setParent(childNode, parent);
    stTree_setBranchLength(childNode, branchLength);
    _nodeMap.insert(pair<string, stTree *>(name, childNode));

    Hdf5Genome *genome = new Hdf5Genome(name, this, _file, _dcprops, _inMemory);
    _openGenomes.insert(pair<string, Hdf5Genome *>(name, genome));
    _dirty = true;
    return genome;
}

Genome *Hdf5Alignment::addRootGenome(const string &name, double branchLength) {
    if (name.empty() == true) {
        throw hal_exception("name can't be empty");
    }
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt != _nodeMap.end()) {
        throw hal_exception(string("node ") + name + " already exists");
    }
    stTree *node = stTree_construct();
    stTree_setLabel(node, name.c_str());
    if (_tree != NULL) {
        stTree_setParent(_tree, node);
        stTree_setBranchLength(_tree, branchLength);
    }
    _tree = node;
    _nodeMap.insert(pair<string, stTree *>(name, node));

    Hdf5Genome *genome = new Hdf5Genome(name, this, _file, _dcprops, _inMemory);
    _openGenomes.insert(pair<string, Hdf5Genome *>(name, genome));
    _dirty = true;
    return genome;
}

// May only make sense to remove a leaf genome
// (so that's what is done here right now)
void Hdf5Alignment::removeGenome(const string &name) {
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt == _nodeMap.end()) {
        throw hal_exception("node " + name + " does not exist");
    }
    stTree *node = findIt->second;
    if (stTree_getChildNumber(node) != 0) {
        throw hal_exception("node " + name + " has a child");
    }
    if (stTree_getParent(node) != NULL) {
        // The parent will have to be updated to fix its bottom segments
        Genome *parentGenome = openGenome(getParentName(name));
        hal_size_t n = parentGenome->getNumBottomSegments();
        vector<string> childNames = getChildNames(parentGenome->getName());
        hal_index_t removedChildIndex = -1;
        for (size_t i = 0; i < childNames.size(); i++) {
            if (childNames[i] == name) {
                removedChildIndex = i;
                break;
            }
        }
        assert(removedChildIndex != -1);
        stTree_setParent(node, NULL);
        // Copy the old bottom segments into memory. Necessary since it isn't
        // possible to update the bottom array in-place
        // Local typedefs
        typedef struct {
            hal_index_t childIndex;
            bool reversed;
        } ChildInfo_t;
        typedef struct {
            size_t start;
            size_t length;
            ChildInfo_t *children;
            hal_index_t topParseIndex;
        } BottomSegment_t;

        BottomSegment_t *bottomSegments = (BottomSegment_t *)malloc(sizeof(BottomSegment_t) * n);
        assert(bottomSegments != NULL);
        BottomSegmentIteratorPtr oldBotIt = parentGenome->getBottomSegmentIterator();
        for (size_t i = 0; (hal_size_t)oldBotIt->getArrayIndex() < n; oldBotIt->toRight(), i++) {
            bottomSegments[i].start = oldBotIt->getStartPosition();
            bottomSegments[i].length = oldBotIt->getLength();
            bottomSegments[i].children = (ChildInfo_t *)malloc(sizeof(ChildInfo_t) * (childNames.size() - 1));
            assert(bottomSegments[i].children != NULL);
            for (hal_index_t oldChild = 0, newChild = 0; oldChild < (hal_index_t)childNames.size(); oldChild++, newChild++) {
                if (oldChild == removedChildIndex) {
                    // Compensate for the increment that will follow
                    newChild--;
                    continue;
                }
                bottomSegments[i].children[newChild].childIndex = oldBotIt->bseg()->getChildIndex(oldChild);
                bottomSegments[i].children[newChild].reversed = oldBotIt->bseg()->getChildReversed(oldChild);
            }
            bottomSegments[i].topParseIndex = oldBotIt->bseg()->getTopParseIndex();
        }
        // Reset the bottom segments. updateBottomDimensions will change the
        // number of children in the bottom segment array to the correct value
        vector<Sequence::UpdateInfo> newBottomDimensions;

        for (SequenceIteratorPtr seqIt = parentGenome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            Sequence::UpdateInfo info(sequence->getName(), sequence->getNumBottomSegments());
            newBottomDimensions.push_back(info);
        }
        parentGenome->updateBottomDimensions(newBottomDimensions);
        ((Hdf5Genome *)parentGenome)->resetBranchCaches();
        // Copy the bottom segments back
        BottomSegmentIteratorPtr botSegIt = parentGenome->getBottomSegmentIterator();
        for (size_t i = 0; (hal_size_t)botSegIt->getArrayIndex() < n; botSegIt->toRight(), i++) {
            botSegIt->setCoordinates(bottomSegments[i].start, bottomSegments[i].length);
            for (hal_index_t child = 0; child < ((hal_index_t)childNames.size()) - 1; child++) {
                botSegIt->bseg()->setChildIndex(child, bottomSegments[i].children[child].childIndex);
                botSegIt->bseg()->setChildReversed(child, bottomSegments[i].children[child].reversed);
            }
            free(bottomSegments[i].children);
            botSegIt->bseg()->setTopParseIndex(bottomSegments[i].topParseIndex);
        }
        free(bottomSegments);
    }
    map<string, Hdf5Genome *>::iterator mapIt = _openGenomes.find(name);
    if (mapIt != _openGenomes.end()) {
        closeGenome(mapIt->second);
    }
    _file->unlink(name);
    _nodeMap.erase(findIt);
    stTree_destruct(node);
    _dirty = true;
}

const Genome *Hdf5Alignment::openGenome(const string &name) const {
    return const_cast<Hdf5Alignment *>(this)->openGenome(name);
}

Genome *Hdf5Alignment::openGenome(const string &name) {
    map<string, Hdf5Genome *>::iterator mapit = _openGenomes.find(name);
    if (mapit != _openGenomes.end()) {
        return mapit->second;
    }
    Hdf5Genome *genome = NULL;
    if (_nodeMap.find(name) != _nodeMap.end()) {
        genome = new Hdf5Genome(name, this, _file, _dcprops, _inMemory);
        _openGenomes.insert(pair<string, Hdf5Genome *>(name, genome));
    }
    return genome;
}

void Hdf5Alignment::closeGenome(const Genome *genome) const {
    string name = genome->getName();
    map<string, Hdf5Genome *>::iterator mapIt = _openGenomes.find(name);
    if (mapIt == _openGenomes.end()) {
        throw hal_exception("Attempt to close non-open genome.  "
                            "Should not even be possible");
    }
    mapIt->second->write();
    delete mapIt->second;
    _openGenomes.erase(mapIt);

    // reset the parent/child genoem cachces (which store genome pointers to
    // the genome we're closing
    if (name != getRootName()) {
        mapIt = _openGenomes.find(getParentName(name));
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
    vector<string> childNames = getChildNames(name);
    for (size_t i = 0; i < childNames.size(); ++i) {
        mapIt = _openGenomes.find(childNames[i]);
        if (mapIt != _openGenomes.end()) {
            mapIt->second->resetBranchCaches();
        }
    }
}

string Hdf5Alignment::getRootName() const {
    if (_tree == NULL) {
        throw hal_exception("Can't get root name of empty tree");
    }
    return stTree_getLabel(_tree);
}

string Hdf5Alignment::getParentName(const string &name) const {
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt == _nodeMap.end()) {
        throw hal_exception(string("node not found: ") + name);
    }
    stTree *node = findIt->second;
    stTree *parent = stTree_getParent(node);
    if (parent == NULL) {
        return "";
    }
    return stTree_getLabel(parent);
}

void Hdf5Alignment::updateBranchLength(const string &parentName, const string &childName, double length) {
    map<string, stTree *>::iterator findIt = _nodeMap.find(childName);
    if (findIt == _nodeMap.end()) {
        throw hal_exception(string("node ") + childName + " not found");
    }
    stTree *node = findIt->second;
    stTree *parent = stTree_getParent(node);
    if (parent == NULL || parentName != stTree_getLabel(parent)) {
        throw hal_exception(string("edge ") + parentName + "--" + childName + " not found");
    }
    stTree_setBranchLength(node, length);
    _dirty = true;
}

double Hdf5Alignment::getBranchLength(const string &parentName, const string &childName) const {
    map<string, stTree *>::iterator findIt = _nodeMap.find(childName);
    if (findIt == _nodeMap.end()) {
        throw hal_exception(string("node ") + childName + " not found");
    }
    stTree *node = findIt->second;
    stTree *parent = stTree_getParent(node);
    if (parent == NULL || parentName != stTree_getLabel(parent)) {
        throw hal_exception(string("edge ") + parentName + "--" + childName + " not found");
    }
    return stTree_getBranchLength(node);
}

vector<string> Hdf5Alignment::getChildNames(const string &name) const {
    map<string, stTree *>::iterator findIt = _nodeMap.find(name);
    if (findIt == _nodeMap.end()) {
        throw hal_exception(string("node ") + name + " not found");
    }
    stTree *node = findIt->second;
    int32_t numChildren = stTree_getChildNumber(node);
    vector<string> childNames(numChildren);
    for (int32_t i = 0; i < numChildren; ++i) {
        childNames[i] = stTree_getLabel(stTree_getChild(node, i));
    }
    return childNames;
}

vector<string> Hdf5Alignment::getLeafNamesBelow(const string &name) const {
    vector<string> leaves;
    vector<string> children;
    deque<string> bfQueue;
    bfQueue.push_front(name);
    while (bfQueue.empty() == false) {
        string &current = bfQueue.back();
        children = getChildNames(current);
        if (children.empty() == true && current != name) {
            leaves.push_back(current);
        }
        for (size_t i = 0; i < children.size(); ++i) {
            bfQueue.push_front(children[i]);
        }
        bfQueue.pop_back();
    }
    return leaves;
}

hal_size_t Hdf5Alignment::getNumGenomes() const {
    if (_tree == NULL) {
        assert(_nodeMap.empty() == true);
        return 0;
    } else {
        return _nodeMap.size();
    }
}

MetaData *Hdf5Alignment::getMetaData() {
    return _metaData;
}

const MetaData *Hdf5Alignment::getMetaData() const {
    return _metaData;
}

string Hdf5Alignment::getNewickTree() const {
    if (_tree == NULL) {
        return "";
    } else {
        char *treeString = stTree_getNewickTreeString(_tree);
        string returnString(treeString);
        free(treeString);
        return returnString;
    }
}

string Hdf5Alignment::getVersion() const {
    try {
        HDF5DisableExceptionPrinting prDisable;
        _file->openGroup(VersionGroupName);
        HDF5MetaData versionMeta(_file, VersionGroupName);
        if (versionMeta.has(VersionGroupName) == false) {
            throw Exception();
        }
        return versionMeta.get(VersionGroupName);
    } catch (Exception &e) {
        // since there was no version tag at the beginning,
        // we retroactively consider it 0.0
        return "0.0";
    }
}

bool Hdf5Alignment::isReadOnly() const {
    return (_flags & (H5F_ACC_RDWR | H5F_ACC_TRUNC)) == 0;
}

void Hdf5Alignment::writeTree() {
    if (_dirty == false)
        return;

    char *treeString = NULL;
    if (_tree != NULL) {
        treeString = stTree_getNewickTreeString(_tree);
    } else {
        treeString = (char *)malloc(sizeof(char));
        treeString[0] = '\0';
    }
    assert(_file != NULL);
    HDF5MetaData treeMeta(_file, TreeGroupName);
    treeMeta.set(TreeGroupName, treeString);
    free(treeString);
}

void Hdf5Alignment::writeVersion() {
    if (_dirty == false)
        return;

    assert(_file != NULL);
    HDF5MetaData versionMeta(_file, VersionGroupName);
    versionMeta.set(VersionGroupName, HAL_VERSION);
}

static void addNodeToMap(stTree *node, map<string, stTree *> &nodeMap) {
    const char *label = stTree_getLabel(node);
    assert(label != NULL);
    string name(label);
    assert(nodeMap.find(name) == nodeMap.end());
    nodeMap.insert(pair<string, stTree *>(name, node));
    int32_t numChildren = stTree_getChildNumber(node);
    for (int32_t i = 0; i < numChildren; ++i) {
        addNodeToMap(stTree_getChild(node, i), nodeMap);
    }
}

void Hdf5Alignment::loadTree() {
    _nodeMap.clear();
    HDF5MetaData treeMeta(_file, TreeGroupName);
    const string &treeString = treeMeta.get(TreeGroupName);
    if (_tree != NULL) {
        stTree_destruct(_tree);
    }
    if (treeString.empty() == true) {
        _tree = stTree_construct();
    } else {
        _tree = stTree_parseNewickString(const_cast<char *>(treeString.c_str()));
        addNodeToMap(_tree, _nodeMap);
    }
}

void Hdf5Alignment::replaceNewickTree(const string &newNewickString) {
    _nodeMap.clear();
    HDF5MetaData treeMeta(_file, TreeGroupName);
    treeMeta.set(TreeGroupName, newNewickString);
    treeMeta.write();
    loadTree();
}
