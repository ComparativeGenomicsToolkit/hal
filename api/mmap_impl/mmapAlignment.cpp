#include "mmapAlignment.h"
#include "halCLParser.h"
#include "mmapGenome.h"

using namespace hal;
using namespace std;

static const int NAME_HASH_GROWTH_FACTOR = 1024; // allow lots of initial space

MMapAlignment::MMapAlignment(const std::string &alignmentPath, unsigned mode, size_t fileSize)
    : _alignmentPath(alignmentPath), _mode(mode), _fileSize(fileSize), _file(NULL), _data(NULL), _genomeNameHash(NULL),
      _tree(NULL) {
    _file = MMapFile::factory(alignmentPath, mode, fileSize);
    if (mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

MMapAlignment::MMapAlignment(const std::string &alignmentPath, unsigned mode, const CLParser *parser)
    : _alignmentPath(alignmentPath), _mode(halDefaultAccessMode(mode)), _fileSize(0), _file(NULL), _data(NULL),
      _genomeNameHash(NULL), _tree(NULL) {
    initializeFromOptions(parser);
    _file = MMapFile::factory(alignmentPath, _mode, _fileSize);
    if (mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

void MMapAlignment::close() {
    // Free the memory used by all open genomes.
    for (auto kv : _openGenomes) {
        delete kv.second;
    }
    // Close the actual file.
    delete _genomeNameHash;
    _genomeNameHash = NULL;
    _file->close();
}

void MMapAlignment::defineOptions(CLParser *parser, unsigned mode) {
    if (mode & CREATE_ACCESS) {
        parser->addOption("mmapFileSize", "mmap HAL file initial size (in gigabytes)", MMAP_DEFAULT_FILE_SIZE_GB);
    } else if (mode & WRITE_ACCESS) {
        parser->addOption("mmapSizeIncrease", "additional space to reserve at end of file (in gigabytes)", 1);
    }
}

/* initialize class from options */
void MMapAlignment::initializeFromOptions(const CLParser *parser) {
    if (_mode & CREATE_ACCESS) {
        _fileSize = GIGABYTE * parser->get<size_t>("mmapFileSize");
    } else if (_mode & WRITE_ACCESS) {
        // TODO: this causes _fileSize's meaning to be far too
        // overloaded: sometimes (CREATE_ACCESS) it is a requested
        // maximum file size, others (WRITE_ACCESS) it is an offset,
        // and still others (READ_ACCESS) it is totally
        // ignored. Probably the cleanest solution is to just create
        // 3 separate factory functions.
        _fileSize = GIGABYTE * parser->get<size_t>("mmapSizeIncrease");
    }
}

void MMapAlignment::create() {
    _file->allocMem(sizeof(MMapAlignmentData), true);
    _data = static_cast<MMapAlignmentData *>(resolveOffset(_file->getRootOffset(), sizeof(MMapAlignmentData)));
    _data->_numGenomes = 0;
    _data->_genomeNameHashOffset = MMAP_NULL_OFFSET;
}

void MMapAlignment::open() {
    _data = static_cast<MMapAlignmentData *>(resolveOffset(_file->getRootOffset(), sizeof(MMapAlignmentData)));
    if (_data->_genomeNameHashOffset != MMAP_NULL_OFFSET) {
        _genomeNameHash = new MMapPerfectHashTable(_file, _data->_genomeNameHashOffset, NAME_HASH_GROWTH_FACTOR);
    }
    loadTree();
}

MMapGenome *MMapAlignmentData::addGenome(MMapAlignment *alignment, const std::string &name) {
    // FIXME: would be nice to allocate extra space and only move when needed.
    size_t newGenomeArraySize = (_numGenomes + 1) * sizeof(MMapGenomeData);
    size_t newGenomeArrayOffset = alignment->allocateNewArray(newGenomeArraySize);
    MMapGenomeData *newGenomeArray = (MMapGenomeData *)alignment->resolveOffset(newGenomeArrayOffset, newGenomeArraySize);
    if (_numGenomes != 0) {
        // Copy over old genome data.
        MMapGenomeData *oldGenomeArray =
            (MMapGenomeData *)alignment->resolveOffset(_genomeArrayOffset, _numGenomes * sizeof(MMapGenomeData));
        memcpy(newGenomeArray, oldGenomeArray, _numGenomes * sizeof(MMapGenomeData));
    }
    _genomeArrayOffset = newGenomeArrayOffset;

    // Update any existing genome's data entries to point to their new location.
    for (auto &name_genome : alignment->_openGenomes) {
        MMapGenome *genome = name_genome.second;
        genome->updateGenomeArrayBasePtr(newGenomeArray);
    }

    MMapGenomeData *data = newGenomeArray + _numGenomes;
    _numGenomes += 1;
    MMapGenome *genome = new MMapGenome(alignment, data, _numGenomes - 1, name);
    return genome;
}

vector<string> MMapAlignmentData::getGenomeNames(MMapAlignment *alignment) {
    vector<string> names;
    MMapGenomeData *genomeArray =
        static_cast<MMapGenomeData *>(alignment->resolveOffset(_genomeArrayOffset, _numGenomes * sizeof(MMapGenomeData)));
    for (int i = 0; i < _numGenomes; i++) {
        names.push_back(genomeArray[i].getName(alignment));
    }
    return names;
}

/* this only adds name to hash function, index must still be set */
void MMapAlignment::addGenomeToNameHash(const MMapGenome *genome, vector<string> &existingNames) {
    const string &name = genome->getName();
    vector<string> newNames;
    newNames.push_back(name);
    assert(find(existingNames.begin(), existingNames.end(), name) == existingNames.end()); // must not exist

    if (_genomeNameHash == NULL) {
        _genomeNameHash = new MMapPerfectHashTable(_file, _data->_genomeNameHashOffset, NAME_HASH_GROWTH_FACTOR);
    }
    _data->_genomeNameHashOffset = _genomeNameHash->addKeys(newNames, existingNames);
    _genomeNameHash->setIndex(name, genome->getArrayIndex());
}

Genome *MMapAlignment::addLeafGenome(const string &name, const string &parentName, double branchLength) {
    _childNames.clear();
    stTree *parentNode = getGenomeNode(parentName);
    stTree *childNode = stTree_construct();
    stTree_setLabel(childNode, name.c_str());
    stTree_setParent(childNode, parentNode);
    stTree_setBranchLength(childNode, branchLength);
    writeTree();
    vector<string> existingNames = _data->getGenomeNames(this);
    MMapGenome *genome = _data->addGenome(this, name);
    addGenomeToNameHash(genome, existingNames);
    _openGenomes[name] = genome;
    return genome;
}

Genome *MMapAlignment::addRootGenome(const string &name, double branchLength) {
    _childNames.clear();
    stTree *newRoot = stTree_construct();
    stTree_setLabel(newRoot, name.c_str());
    if (_tree != NULL) {
        stTree_setParent(_tree, newRoot);
        stTree_setBranchLength(_tree, branchLength);
    }
    _tree = newRoot;
    writeTree();
    vector<string> existingNames = _data->getGenomeNames(this);
    MMapGenome *genome = _data->addGenome(this, name);
    addGenomeToNameHash(genome, existingNames);
    _openGenomes[name] = genome;
    return genome;
}

Genome *MMapAlignment::_openGenome(const string &name) const {
    if (_openGenomes.find(name) != _openGenomes.end()) {
        // Already loaded.
        return _openGenomes[name];
    }
    if (_genomeNameHash == NULL) {
        throw GenomeNotFoundException(name);
    }
    hal_index_t genomeIndex = _genomeNameHash->getIndex(name);
    if (genomeIndex == NULL_INDEX) {
        throw GenomeNotFoundException(name);
    }
    // FIXME: put this in a function
    MMapGenomeData *genomeDataArray =
        (MMapGenomeData *)resolveOffset(_data->_genomeArrayOffset, _data->_numGenomes * sizeof(MMapGenomeData));
    if (genomeDataArray[genomeIndex].getName(const_cast<MMapAlignment *>(this)) != name) {
        throw GenomeNotFoundException(name); // name not in perfect hash
    }
    MMapGenome *genome = new MMapGenome(const_cast<MMapAlignment *>(this), &genomeDataArray[genomeIndex], genomeIndex);
    _openGenomes[name] = genome;
    return genome;
}
