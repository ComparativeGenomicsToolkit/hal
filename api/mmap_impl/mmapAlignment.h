#ifndef _MMAPALIGNMENT_H
#define _MMAPALIGNMENT_H
#include <deque>
#include "halAlignment.h"
#include "mmapFile.h"
#include "sonLib.h"

namespace hal {
class MMapAlignment;
class MMapGenome;
class MMapAlignmentData {
    friend class MMapAlignment;
public:
    const char *getNewickString(const MMapAlignment *alignment);
    void setNewickString(MMapAlignment *alignment, const char *newickString);

    MMapGenome *addGenome(MMapAlignment *alignment, const std::string &name);
protected:
    size_t _numGenomes;
private:
    size_t _newickStringOffset;
    size_t _newickStringLength;
    size_t _genomeArrayOffset;
};

class MMapAlignment : public Alignment {
    public:
    /* constructor with all arguments specified */
    MMapAlignment(const std::string& alignmentPath,
                  unsigned mode = READ_ACCESS,
                  size_t initSize = MMAP_DEFAULT_INIT_SIZE,
                  size_t growSize = MMAP_DEFAULT_GROW_SIZE);
    
    /* constructor from command line options */
    MMapAlignment(const std::string& alignmentPath,
                  unsigned mode,
                  CLParserConstPtr parser);

    // Allocate new array and return the offset.
    size_t allocateNewArray(size_t size) const { return _file->allocMem(size, false); };
    void *resolveOffset(size_t offset, size_t len) const { return _file->toPtr(offset, len); };

    void close() { _file->close(); }

    const std::string& getStorageFormat() const {
        return STORAGE_FORMAT_MMAP;
    }

    Genome* addLeafGenome(const std::string& name,
                          const std::string& parentName,
                          double branchLength);

    Genome* addRootGenome(const std::string& name,
                          double branchLength);

    void removeGenome(const std::string& name) {
        // FIXME
        throw hal_exception("unimplemented; don't want to deal with this right now.");
    };

    Genome* insertGenome(const std::string& name,
                         const std::string& parentName,
                         const std::string& childName,
                         double upperBranchLength) {
        // FIXME
        throw hal_exception("unimplemented; don't want to deal with this right now.");
    };

    const Genome* openGenome(const std::string& name) const {
        Genome *genome = _openGenome(name);
        return const_cast<const Genome *>(genome);
    };

    Genome* openGenome(const std::string& name) {
        return _openGenome(name);
    }

    void closeGenome(const Genome* genome) const {
        // Intentionally left blank.
    };

    std::string getRootName() const {
        return stTree_getLabel(_tree);
    };

    std::string getParentName(const std::string& name) const {
        stTree *node = getGenomeNode(name);
        stTree *parent = stTree_getParent(node);
        if (parent == NULL) {
            throw hal_exception("Root genome has no parent");
        }
        return stTree_getLabel(parent);
    };

    void updateBranchLength(const std::string& parentName,
                            const std::string& childName,
                            double length) {
        stTree *node = getGenomeNode(childName);
        stTree *parent = stTree_getParent(node);
        if (stTree_getLabel(parent) != parentName) {
            throw hal_exception("Branch doesn't exist.");
        }
        stTree_setBranchLength(node, length);
        writeTree();
    }

    double getBranchLength(const std::string& parentName,
                           const std::string& childName) const {
        stTree *node = getGenomeNode(childName);
        stTree *parent = stTree_getParent(node);
        if (stTree_getLabel(parent) != parentName) {
            throw hal_exception("Branch doesn't exist.");
        }
        return stTree_getBranchLength(node);
    };

    std::vector<std::string> 
    getChildNames(const std::string& name) const {
        stTree *node = getGenomeNode(name);
        std::vector<std::string> ret;
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            const char *name = stTree_getLabel(stTree_getChild(node, i));
            ret.push_back(std::string(name));
        }
        return ret;
    };

    std::vector<std::string>
    getLeafNamesBelow(const std::string& name) const {
        std::vector<std::string> leaves;
        std::vector<std::string> children;
        std::deque<std::string> bfQueue;
        bfQueue.push_front(name);
        while (bfQueue.empty() == false)
        {
            std::string& current = bfQueue.back();
            children = getChildNames(current);
            if (children.empty() == true && current != name)
            {
                leaves.push_back(current);
            }
            for (size_t i = 0; i < children.size(); ++i)
            {
                bfQueue.push_front(children[i]);
            }
            bfQueue.pop_back();
        }
        return leaves;
    };

    hal_size_t getNumGenomes() const { return _data->_numGenomes; };

    MetaData* getMetaData() {
        throw hal_exception("unimplemented");
    };

    const MetaData* getMetaData() const {
        throw hal_exception("unimplemented");
    };

    std::string getNewickTree() const { return _data->getNewickString(this); };

    std::string getVersion() const { return _file->getVersion(); };

    bool isReadOnly() const { return _file->isReadOnly(); };

    void replaceNewickTree(const std::string &newNewickString) {
        _data->setNewickString(this, newNewickString.c_str());
        loadTree();
    };

private:
    void create();
    void open();
    Genome *_openGenome(const std::string &name) const;
    stTree *getGenomeNode(const std::string &name) const {
        stTree *node = stTree_findChild(_tree, name.c_str());
        if (node == NULL) {
            throw hal_exception("genome " + name + " not found in alignment.");
        }
        return node;
    };
    void loadTree() {
        if (_tree != NULL) {
            stTree_destruct(_tree);
        }
        _tree = stTree_parseNewickString(_data->getNewickString(this));
    };
    void writeTree() {
        char *newickString = stTree_getNewickTreeString(_tree);
        _data->setNewickString(this, newickString);
        free(newickString);
    };
    MMapFile *_file;
    MMapAlignmentData *_data;
    stTree *_tree;
    mutable std::map<std::string, MMapGenome *> _openGenomes;
};

inline const char *MMapAlignmentData::getNewickString(const MMapAlignment *alignment) {
    return (const char *) alignment->resolveOffset(_newickStringOffset, _newickStringLength);
}

inline void MMapAlignmentData::setNewickString(MMapAlignment *alignment, const char *newickString) {
    _newickStringLength = strlen(newickString) + 1;
    _newickStringOffset = alignment->allocateNewArray(_newickStringLength);
    char *dest = (char *) alignment->resolveOffset(_newickStringOffset, _newickStringLength);
    strncpy(dest, newickString, _newickStringLength);
}
}
#endif
