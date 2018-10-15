#include "mmapAlignment.h"
#include "mmapGenome.h"
#include "mmapCLParser.h"

using namespace hal;
using namespace std;


MMapAlignment::MMapAlignment(const std::string& alignmentPath,
                             unsigned mode,
                             size_t initSize,
                             size_t growSize):
    _file(NULL), _data(NULL), _tree(NULL) {
    _file = MMapFile::localFactory(alignmentPath, mode, initSize, growSize);
    if (mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
}

MMapAlignment::MMapAlignment(const std::string& alignmentPath,
                             unsigned mode,
                             CLParserConstPtr parser):
    _file(NULL), _data(NULL), _tree(NULL) {
    _file = MMapFile::localFactory(alignmentPath, mode,
                                   MMapCLParser::getInitSize(parser),
                                   MMapCLParser::getGrowSize(parser));
    if (mode & CREATE_ACCESS) {
        create();
    } else {
        open();
    }
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
    MMapGenomeData *data = newGenomeArray + _numGenomes;
    _numGenomes += 1;
    MMapGenome *genome = new MMapGenome(alignment, data, name);
    return genome;
}

Genome* MMapAlignment::addLeafGenome(const string& name,
                                     const string& parentName,
                                     double branchLength) {
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
        MMapGenome curGenome(const_cast<MMapAlignment *>(this), genomeArray + i);
        if (curGenome.getName() == name) {
            genome = new MMapGenome(const_cast<MMapAlignment *>(this), genomeArray + i);
            break;
        }
    }
    if (genome == NULL) {
        throw hal_exception("Genome " + name + "not found in alignment");
    }
    _openGenomes[name] = genome;
    return genome;
}
