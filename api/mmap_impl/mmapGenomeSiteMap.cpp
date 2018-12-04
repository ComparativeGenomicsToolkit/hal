#include "mmapGenomeSiteMap.h"
#include "mmapRbTree.h"
#include "mmapSequence.h"
using namespace std;
using namespace hal;

namespace hal {
    /* Object to store in tree, key is pointer to MMapGenomeSiteMapNode.
     * Since the rbtree library is no allocation, this also holds the nodes.
     * We initialize a vector of these at the start. */
    class TmpTreeNode {
        public:
        MMapGenomeSiteMapNode _siteMapNode;
        struct rb_tree_node _rbTreeNode;
    };

    class TmpTreeNodes: public vector<TmpTreeNode> {
    };
    
    
    /* comparison function for MMapGenomeSiteMapNode  */
    static int mmapGenomeSiteMapNodeCmp(const void *lhs, const void *rhs) {
        const MMapGenomeSiteMapNode* lnode = static_cast<const MMapGenomeSiteMapNode*>(lhs);
        const MMapGenomeSiteMapNode* rnode = static_cast<const MMapGenomeSiteMapNode*>(rhs);
        if ((lnode->_startPosition + lnode->_length) <= rnode->_startPosition) {
            return -1;
        } else if (lnode->_startPosition >= (rnode->_startPosition + rnode->_length)) {
            return 1;
        } else {
            return 0;
        }
    }
}


/* calculate space required for the map in bytes */
size_t hal::MMapGenomeSiteMap::calcRequiredSpace(size_t numSequences) {
    // root is in header
    return MMapFile::alignRound(sizeof(MMapGenomeSiteMapData))
        + ((numSequences - 1) * MMapFile::alignRound(sizeof(MMapGenomeSiteMapNode)));
}

/* read header information */
void hal::MMapGenomeSiteMap::readGsm(size_t gsmOffset) {
    _gsmOffset = gsmOffset;
    _data = static_cast<MMapGenomeSiteMapData*>(_file->toPtr(gsmOffset, sizeof(MMapGenomeSiteMapData)));
    // prefetch full table
    _file->toPtr(gsmOffset, calcRequiredSpace(_data->_numSequences));
}

void hal::MMapGenomeSiteMap::createGsm(size_t numSequences) {
    _gsmOffset = _file->allocMem(calcRequiredSpace(numSequences));
    _data = static_cast<MMapGenomeSiteMapData*>(_file->toPtr(_gsmOffset, calcRequiredSpace(numSequences)));
    _data->_numSequences = numSequences;
}

/* load into temporary site tree, which is balanced */
void hal::MMapGenomeSiteMap::loadTmpTree(const vector<MMapSequence*>& sequences,
                                         struct rb_tree *tmpTree,
                                         TmpTreeNodes& tmpTreeNodes) {
    tmpTreeNodes.resize(sequences.size());  // preallocate so memory doesn't change
    int iNode = 0;
    for (auto seq: sequences) {
        TmpTreeNode* node = &tmpTreeNodes[iNode++];
        node->_siteMapNode = MMapGenomeSiteMapNode(seq->getStartPosition(),
                                                   seq->getSequenceLength(),
                                                   seq->getArrayIndex());
        rb_result_t stat = rb_tree_insert(tmpTree, &node->_siteMapNode, &node->_rbTreeNode);
        if (stat != RB_OK) {
            throw hal_exception("error inserting in rb_tree");
        }
    }
}

/* recursively copy tmp tree to mmap file */
hal_index_t hal::MMapGenomeSiteMap::copyTree(struct rb_tree_node *tmpNode,
                                             int& nextNodeIdx) {
    if (tmpNode == NULL) {
        return NULL_INDEX;
    }
    const MMapGenomeSiteMapNode *tmpGsmNode = reinterpret_cast<const MMapGenomeSiteMapNode *>(tmpNode->key);
    // allocate and copy
    hal_index_t nodeIdx = nextNodeIdx++;
    MMapGenomeSiteMapNode* gsmNode = getNodePtr(nodeIdx);
    *gsmNode = *tmpGsmNode;
    
    // children
    gsmNode->_leftNodeIndex = copyTree(tmpNode->left, nextNodeIdx);
    gsmNode->_rightNodeIndex = copyTree(tmpNode->right, nextNodeIdx);
    
    return nodeIdx;
}

size_t hal::MMapGenomeSiteMap::build(const vector<MMapSequence*>& sequences) {
    struct rb_tree tmpTree;
    TmpTreeNodes tmpTreeNodes;  // manages memory for tmp tree

    rb_tree_new(&tmpTree, mmapGenomeSiteMapNodeCmp);
    loadTmpTree(sequences, &tmpTree, tmpTreeNodes);
    createGsm(sequences.size());
    int nextNodeIdx = 0;
    copyTree(tmpTree.root, nextNodeIdx);
    return _gsmOffset;
}


hal_index_t MMapGenomeSiteMap::getSequenceIndexBySite(size_t position) {
    assert(_gsmOffset != MMAP_NULL_OFFSET);
    const MMapGenomeSiteMapNode* node = getNodePtr(0);
    while (node != NULL) {
        int dir = node->positionCmp(position);
        if (dir < 0) {
            node = getNodePtr(node->_leftNodeIndex);
        } else if (dir > 0) {
            node = getNodePtr(node->_rightNodeIndex);
        } else {
            return node->_sequenceIndex;
        }
    }
    return NULL_INDEX;
}
