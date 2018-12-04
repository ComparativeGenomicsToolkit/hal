#ifndef _MMAPGENOMESITEMAP_h
#define _MMAPGENOMESITEMAP_h
#include "mmapFile.h"
#include <string>
#include <vector>
struct rb_tree;
struct rb_tree_node;

namespace hal {
    class MMapSequence;
    class TmpTreeNodes;
    
    /* node in the binary search tree */
    class MMapGenomeSiteMapNode {
        public:
        MMapGenomeSiteMapNode(size_t startPosition = 0,
                              size_t length = 0,
                              hal_index_t sequenceIndex = NULL_INDEX):
            _startPosition(startPosition),
            _length(length),
            _sequenceIndex(sequenceIndex),
            _leftNodeIndex(NULL_INDEX),
            _rightNodeIndex(NULL_INDEX) {
        }

        /* check if positions is in entry?  Return 0 if contain, -1 if less than, or
         * 1 if great than */
        int positionCmp(size_t position) const {
            if (position < _startPosition) {
                return -1;
            } else if (position >= (_startPosition + _length)) {
                return 1;
            } else {
                return 0;
            }
        }
        
        
        /* sequence range in genome */
        size_t _startPosition;
        size_t _length;
        hal_index_t _sequenceIndex;

        /* child links */
        hal_index_t _leftNodeIndex;
        hal_index_t _rightNodeIndex;
    };
    
    class MMapGenomeSiteMapData {
        public:
        size_t _numSequences;
        // other nodes immediately follow
        MMapGenomeSiteMapNode _root; 
    };
    
    /**
     * MMap file structure used to map position in genome to specific
     * sequence.  This builds a balance binary tree and stores it in the
     * mmapped file for direct access.
     */
    class MMapGenomeSiteMap {
        public:
        /** Construct new object for accessing site map in HAL file.
         * If the hash table is being created, then gsmOffset
         * should be MMAP_NULL_OFFSET.  */
        MMapGenomeSiteMap(MMapFile* mmapFile,
                          size_t gsmOffset):
            _file(mmapFile),
            _gsmOffset(gsmOffset),
            _data(NULL) {
            if (gsmOffset != MMAP_NULL_OFFSET) {
                readGsm(gsmOffset);
            }
        }

        /* Build the map from an array of sequences, return gsmOffset.
         * Rebuilding will just lose space in the file. */
        size_t build(const std::vector<MMapSequence*>& sequences);

        /** find the sequence index containing a position */
        hal_index_t getSequenceIndexBySite(size_t position);

        /** find the sequence containing a position */
        const Sequence* getSequenceBySite(hal_size_t position) const {
            return const_cast<MMapGenomeSiteMap*>(this)->getSequenceBySite(position);
        }

        
        private:
        static size_t calcRequiredSpace(size_t numSequences);
        void readGsm(size_t gsmOffset);
        void createGsm(size_t numSequences);
        void loadTmpTree(const std::vector<MMapSequence*>& sequences,
                         struct rb_tree *tmpTree,
                         TmpTreeNodes& tmpTreeNodes);
        hal_index_t copyTree(struct rb_tree_node *tmpNode,
                             int& nextNodeIdx);

        /* returns null for NULL_INDEX */
        MMapGenomeSiteMapNode* getNodePtr(int nodeIndex) {
            if (nodeIndex == NULL_INDEX) {
                return NULL;
            } else {
                return &_data->_root + nodeIndex;
            }
        }

        /* returns null for NULL_INDEX */
        const MMapGenomeSiteMapNode* getNodePtr(int nodeIndex) const {
            return const_cast<MMapGenomeSiteMap*>(this)->getNodePtr(nodeIndex);
        }
        
        MMapFile* _file;
        size_t _gsmOffset;
        MMapGenomeSiteMapData *_data;
    };
}
#endif
