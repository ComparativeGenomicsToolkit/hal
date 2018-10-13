#include "mmapAlignment.h"
#include "mmapGenome.h"

using namespace hal;
using namespace std;

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
    _numGenomes += 1;
    MMapGenomeData *data = newGenomeArray + _numGenomes;
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
    if (_tree != NULL) {
        // FIXME
        throw hal_exception("adding a new root when one is already defined is unimplemented");
    } else {
        _tree = stTree_construct();
        stTree_setLabel(_tree, name.c_str());
        stTree_setBranchLength(_tree, branchLength);
        MMapGenome *genome = _data->addGenome(this, name);
        _openGenomes[name] = genome;
        return genome;
    }
}

Genome *MMapAlignment::_openGenome(const string &name) const {
    if (_openGenomes.find(name) != _openGenomes.end()) {
        // Already loaded.
        return _openGenomes[name];
    }
    // TODO go through genome array checking names (slow, but necessary for now)
}
