#include "mmapAlignment.h"
#include "mmapGenome.h"
#include "halCLParser.h"

using namespace hal;
using namespace std;


MMapAlignment::MMapAlignment(const std::string& alignmentPath,
                             unsigned mode,
                             size_t initSize,
                             size_t growSize):
    _alignmentPath(alignmentPath),
    _mode(mode), _initSize(initSize), _growSize(initSize),
    _file(NULL), _data(NULL), _tree(NULL) {
    _file = MMapFile::factory(alignmentPath, mode, initSize, growSize);
    if (mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

MMapAlignment::MMapAlignment(const std::string& alignmentPath,
                             unsigned mode,
                             const CLParser* parser):
    _alignmentPath(alignmentPath),
    _mode(halDefaultAccessMode(mode)), _initSize(0), _growSize(0),
    _file(NULL), _data(NULL), _tree(NULL) {
    initializeFromOptions(parser);
    _file = MMapFile::factory(alignmentPath, _mode, _initSize, _growSize, _udcCacheDir);
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
    _file->close();
}

void MMapAlignment::defineOptions(CLParser* parser,
                                  unsigned mode) {
    if (mode & CREATE_ACCESS) {
      parser->addOption("mmapInitSize", "mmap HAL file initial size", MMAP_DEFAULT_INIT_SIZE);
    }
    if (mode & (CREATE_ACCESS | WRITE_ACCESS)) {
      parser->addOption("mmapGrowSize", "mmap HAL file size to grow when more memory is needed", MMAP_DEFAULT_INIT_SIZE);
    }
}

/* initialize class from options */
void MMapAlignment::initializeFromOptions(const CLParser* parser) {
    if (_mode & CREATE_ACCESS) {
        _initSize = parser->get<size_t>("mmapInitSize");
    }
    if (_mode & WRITE_ACCESS) {
        _growSize = parser->get<size_t>("mmapGrowSize");
    }
#ifdef ENABLE_UDC
    if ((_mode & WRITE_ACCESS) == 0) {
        _udcCacheDir = parser->getOption<const string&>("udcCacheDir");
    }
#endif
}

void MMapAlignment::create() {
    _file->allocMem(sizeof(MMapAlignmentData), true);
    _data = static_cast<MMapAlignmentData *>(resolveOffset(_file->getRootOffset(), sizeof(MMapAlignmentData)));
    _data->_numGenomes = 0;
}

void MMapAlignment::open() {
    _data = static_cast<MMapAlignmentData *>(resolveOffset(_file->getRootOffset(), sizeof(MMapAlignmentData)));
    loadTree();
}


MMapGenome *MMapAlignmentData::addGenome(MMapAlignment *alignment, const std::string &name) {
    size_t newGenomeArraySize = (_numGenomes + 1) * sizeof(MMapGenomeData);
    size_t newGenomeArrayOffset = alignment->allocateNewArray(newGenomeArraySize);
    MMapGenomeData *newGenomeArray = (MMapGenomeData *) alignment->resolveOffset(newGenomeArrayOffset, newGenomeArraySize);
    if (_numGenomes != 0) {
        // Copy over old genome data.
        MMapGenomeData *oldGenomeArray = (MMapGenomeData *) alignment->resolveOffset(_genomeArrayOffset, _numGenomes * sizeof(MMapGenomeData));
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

Genome* MMapAlignment::addLeafGenome(const string& name,
                                     const string& parentName,
                                     double branchLength) {
    _childNames.clear();
    stTree *parentNode = getGenomeNode(parentName);
    stTree *childNode = stTree_construct();
    stTree_setLabel(childNode, name.c_str());
    stTree_setParent(childNode, parentNode);
    stTree_setBranchLength(childNode, branchLength);
    writeTree();
    MMapGenome *genome = _data->addGenome(this, name);
    _openGenomes[name] = genome;
    return genome;
}

Genome* MMapAlignment::addRootGenome(const string& name,
                                     double branchLength) {
    _childNames.clear();
    stTree *newRoot = stTree_construct();
    stTree_setLabel(newRoot, name.c_str());
    if (_tree != NULL) {
        stTree_setParent(_tree, newRoot);
        stTree_setBranchLength(_tree, branchLength);
    }
    _tree = newRoot;
    writeTree();
    MMapGenome *genome = _data->addGenome(this, name);
    _openGenomes[name] = genome;
    return genome;
}

Genome *MMapAlignment::_openGenome(const string &name) const {
    if (_openGenomes.find(name) != _openGenomes.end()) {
        // Already loaded.
        return _openGenomes[name];
    }
    // Go through the genome array and find the genome we want.
    MMapGenomeData *genomeArray = (MMapGenomeData *) resolveOffset(_data->_genomeArrayOffset, _data->_numGenomes * sizeof(MMapGenomeData));
    MMapGenome *genome = NULL;
    for (size_t i = 0; i < _data->_numGenomes; i++) {
        MMapGenome curGenome(const_cast<MMapAlignment *>(this), genomeArray + i, i);
        if (curGenome.getName() == name) {
            genome = new MMapGenome(const_cast<MMapAlignment *>(this), genomeArray + i, i);
            break;
        }
    }
    if (genome != NULL) {
        // This genome exists. Save it so we only "open" it once.
        _openGenomes[name] = genome;
    }
    return genome;
}
