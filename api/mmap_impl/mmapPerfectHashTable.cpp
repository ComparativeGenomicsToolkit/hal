#include "mmapPerfectHashTable.h"
#include <algorithm>
#include <cerrno>
#include <memory>
#include <string.h>

using namespace std;
using namespace hal;

/* Fixed values to use for hash function generation */
static const size_t DEFAULT_PHF_LAMBDA = 1; // keys per bucket
static const size_t DEFAULT_PHF_ALPHA = 80; // load factor
static const phf_seed_t DEFAULT_PHF_SEED = 0;
static const bool DEFAULT_PHF_NODIV = true;

/* get relative offset of displacement map */
size_t hal::MMapPerfectHashTable::displacementMapRelOffset() const {
    return MMapFile::alignRound(sizeof(PerfectHashTableData));
}

/* get relative offset of hash table */
size_t hal::MMapPerfectHashTable::hashTableRelOffset(const struct phf *phf) const {
    return MMapFile::alignRound(sizeof(PerfectHashTableData)) + MMapFile::alignRound(phf->r * sizeof(uint32_t));
}

/* compute require space for storing hash given contents of _phf */
size_t hal::MMapPerfectHashTable::calcRequiredSpace(const struct phf *phf) const {
    return MMapFile::alignRound(sizeof(PerfectHashTableData)) + MMapFile::alignRound(phf->r * sizeof(uint32_t)) +
           MMapFile::alignRound(phf->m * sizeof(size_t));
}

/* initialize phf object and this object from file */
void hal::MMapPerfectHashTable::readPhf(size_t phtOffset) {
    // FIXME: change this to do only two prefetches by specifying full length
    _phtOffset = phtOffset;
    _data = static_cast<PerfectHashTableData *>(_file->toPtr(phtOffset, _file->alignRound(sizeof(PerfectHashTableData))));
    _phf.nodiv = _data->_phf_nodiv;
    _phf.seed = _data->_phf_seed;
    _phf.r = _data->_phf_r;
    _phf.m = _data->_phf_m;
    _phf.g = static_cast<uint32_t *>(_file->toPtr(_data->_gOffset, _data->_phf_r * sizeof(uint32_t)));
    _phf.d_max = _data->_phf_d_max;
    _phf.g_op = static_cast<enum phf::g_op_t>(_data->_phf_g_op);
    _hashTable = static_cast<hal_index_t *>(_file->toPtr(_data->_hashTableOffset, _data->_phf_m * sizeof(size_t)));
}

/* Write a phf object fields to file */
void hal::MMapPerfectHashTable::writePhf(size_t phtOffset, struct phf *phf) {
    _data->_phf_nodiv = phf->nodiv;
    _data->_phf_seed = phf->seed;
    _data->_phf_r = phf->r;
    _data->_phf_m = phf->m;
    _data->_phf_d_max = phf->d_max;
    _data->_phf_g_op = phf->g_op;
    _data->_gOffset = phtOffset + displacementMapRelOffset();
    _data->_hashTableOffset = phtOffset + hashTableRelOffset(phf);

    // Calculate the size of the allocated region (the element size is
    // different under different conditions). Keep in sync with
    // PHF::compact().
    size_t size = _data->_phf_r;
    if (_data->_phf_d_max <= 255) {
        size *= sizeof(uint8_t);
    } else if (_data->_phf_d_max <= 65535) {
        size *= sizeof(uint16_t);
    } else {
        size *= sizeof(uint32_t);
    }

    uint32_t *g = static_cast<uint32_t *>(_file->toPtr(_data->_gOffset, size));
    memcpy(g, phf->g, size);

    // initialize hash table to NULL_INDEX
    _hashTable = static_cast<hal_index_t *>(_file->toPtr(_data->_hashTableOffset, _data->_phf_m * sizeof(size_t)));
    uninitialized_fill(_hashTable, _hashTable + _data->_phf_m, NULL_INDEX);
}

vector<hal_index_t> hal::MMapPerfectHashTable::readExistingIndexes(const vector<string> &existingKeys) const {
    vector<hal_index_t> existingIndexes;
    for (auto const &key : existingKeys) {
        existingIndexes.push_back(getIndex(key));
    }
    return existingIndexes;
}

void hal::MMapPerfectHashTable::writeExistingIndexes(const vector<string> &existingKeys,
                                                     const vector<hal_index_t> &existingIndexes) {
    for (int i = 0; i < existingKeys.size(); i++) {
        setIndex(existingKeys[i], existingIndexes[i]);
    }
}

/* Initial or reinitialize space in HAL file. Allocate space if needed,
possibly relocating */
void hal::MMapPerfectHashTable::setupSpace(const struct phf *newPhf) {
    size_t currentUsedSpace = 0;
    size_t currentAllocatedSpace = 0;
    if (_data != NULL) {
        currentUsedSpace = calcRequiredSpace(&_phf);
        currentAllocatedSpace = _data->_allocatedSpace;
    }

    // clear current space, even if moving, so unused is always zero.
    if (_data != NULL) {
        memset(_data, 0, currentUsedSpace);
    }

    // allocate new space if needed
    size_t newNeededSpace = calcRequiredSpace(newPhf);
    if ((_data == NULL) || (newNeededSpace > currentAllocatedSpace)) {
        size_t newAllocatedSpace = newNeededSpace + (2 * _growthFactor * newNeededSpace);
        _phtOffset = _file->allocMem(newAllocatedSpace);
        _data = static_cast<hal::PerfectHashTableData *>(_file->toPtr(_phtOffset, newAllocatedSpace));
        _data->_allocatedSpace = newAllocatedSpace;
    } else {
        _data->_allocatedSpace = currentAllocatedSpace;
    }

    _data->_gOffset = displacementMapRelOffset();
    _data->_hashTableOffset = hashTableRelOffset(newPhf);
}

void hal::MMapPerfectHashTable::buildHash(const vector<string> &newKeys, const vector<string> &existingKeys) {
    vector<string> allKeys;
    std::move(existingKeys.begin(), existingKeys.end(), std::back_inserter(allKeys));
    std::move(newKeys.begin(), newKeys.end(), std::back_inserter(allKeys));

    // build new function in heap
    struct phf newPhf;
    memset(&newPhf, 0, sizeof(newPhf));
    phf_error_t err = PHF::init<string, DEFAULT_PHF_NODIV>(&newPhf, &allKeys[0], allKeys.size(), DEFAULT_PHF_LAMBDA,
                                                           DEFAULT_PHF_ALPHA, DEFAULT_PHF_SEED);
    if (err != 0) {
        throw hal_exception("can't create perfect hash function: " + string(strerror(err)));
    }
    PHF::compact(&newPhf);

    // store in file
    setupSpace(&newPhf);
    writePhf(_phtOffset, &newPhf);

    // reload cached fields
    readPhf(_phtOffset);
    PHF::destroy(&newPhf);
}

size_t hal::MMapPerfectHashTable::addKeys(const vector<string> &newKeys, const vector<string> &existingKeys) {
    vector<hal_index_t> existingIndexes = readExistingIndexes(existingKeys);
    buildHash(newKeys, existingKeys);
    writeExistingIndexes(existingKeys, existingIndexes);
    return _phtOffset;
}
