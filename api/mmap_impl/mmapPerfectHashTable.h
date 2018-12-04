#ifndef _MMAPPERFECTHASHTABLE_H
#define _MMAPPERFECTHASHTABLE_H
#include <string>
#include <vector>
#include "halDefs.h"
#include "mmapFile.h"
#include "mmapPhf.h"

using namespace std;

namespace hal {
    class MMapFile;

    class PerfectHashTableData {
    public:
        /* values from phf object stored as-is. */
        bool _phf_nodiv;
        phf_seed_t _phf_seed;
	size_t _phf_r;           // number of elements in g
	size_t _phf_m;           // number of elements in perfect hash
	size_t _phf_d_max;        // maximum displacement value in g
        uint32_t _phf_g_op;       // operator

        /* transformed pointers to offset */
	size_t _gOffset;           // offset in file of displacement map
        size_t _hashTableOffset;   // offset in file of hash table

        /* extra fields */
        size_t _allocatedSpace;     // total amount of space allocated at this location
    };

    /**
     * Minimal perfect hash table of unique strings stored in an mmap format
     * HAL file to hal_index_t.  Incremental additions are supported by
     * allocating extra memory and rehashing the old sequences with the new
     * one.  If there is insufficient memory, the table is moved.
     */
    class MMapPerfectHashTable {
        public:
        /** Construct new object for accessing hash table in HAL file.
         * If the hash table is being created, then phtOffset
         * should be MMAP_NULL_OFFSET.  In this cases, the table is create
         * when keys are first added.  The growthFactor specified the approximate
         * number of additional set of keys will be added to the file.  Specifies
         * 0 if all keys are added in one set of keys.  This allows for space to
         * avoid moving the hash. */
        MMapPerfectHashTable(MMapFile* mmapFile,
                             size_t phtOffset,
                             size_t growthFactor = 0) :
            _file(mmapFile),
            _phtOffset(phtOffset),
            _data(NULL),
            _growthFactor(growthFactor),
            _hashTable(NULL) {
            memset(&_phf, 0, sizeof(_phf));
            if (phtOffset != MMAP_NULL_OFFSET) {
                readPhf(phtOffset);
            }
        }

        /** Given a string, load the corresponding hal_index.  If key is not in
         * the hash, then either NULL_INDEX or some other offset is
         * returned.  It is up to the caller to validated that the key
         * matches.
         */
        hal_index_t getIndex(const string& key) const {
            return _hashTable[PHF::hash(&_phf, key)];
        }

        /** Save mapping of string o hal_index_t in the hash.  If the string
         * is not in the hash, bad things will happen without notice */
        void setIndex(const string& key,
                      hal_index_t index) {
            _hashTable[PHF::hash(&_phf, key)] = index;
        }

        /** Add a set of keys.  Since all keys must be defined up front, if
         * there are existing keys, then these must be resupplied.  The offset
         * to the table is returned.  If the table new or incrementally change, this
         * maybe different than the original offset and should be saved. */
        size_t addKeys(const vector<string>& newKeys,
                       const vector<string>& existingKeys);

        /** Add a set of keys without any existing keys */
        size_t addKeys(const vector<string>& newKeys) {
            vector<string> existingKeys;
            return addKeys(newKeys, existingKeys);
        }
        
        private:
        inline size_t displacementMapRelOffset() const;
        inline size_t hashTableRelOffset(const struct phf* phf) const;
        inline size_t calcRequiredSpace(const struct phf* phf) const;
        void readPhf(size_t phtOffset);
        void writePhf(size_t phtOffset,
                      struct phf* phf);
        vector<hal_index_t> readExistingIndexes(const vector<string>& existingKeys) const;
        void writeExistingIndexes(const vector<string>& existingKeys,
                                  const vector<hal_index_t>& existingIndexes);
        void setupSpace(const struct phf *newPhf);
        void buildHash(const vector<string>& newKeys,
                       const vector<string>& existingKeys);
        
        MMapFile* _file;
        size_t _phtOffset;  // pht is perfect hash table, to distinguish from phf
        PerfectHashTableData *_data;
        size_t _growthFactor;
        struct phf _phf;
        hal_index_t *_hashTable;
    };
}

#endif

// Local Variables:
// mode: c++
// End:
