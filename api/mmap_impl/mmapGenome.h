#ifndef _MMAPGENOME_H
#define _MMAPGENOME_H
#include <map>
#include "halGenome.h"
#include "mmapAlignment.h"

namespace hal {
class MMapTopSegmentData;
class MMapBottomSegmentData;
class MMapSequence;

class MMapGenomeData {
    friend class MMapGenome;
protected:
    hal_size_t _nameLength;
    hal_size_t _totalSequenceLength;
    hal_size_t _numSequences;
    hal_size_t _numMetadata;
    hal_size_t _numTopSegments;
    hal_size_t _numBottomSegments;
    hal_size_t _numChildren;

    // All offsets are relative to the start of this MMapGenomeData, in
    // bytes. Negative offsets are allowed, but unlikely.
    int64_t _nameOffset;
    int64_t _sequencesOffset;
    int64_t _metadataOffset;
    int64_t _dnaOffset;
    int64_t _topSegmentsOffset;
    int64_t _bottomSegmentsOffset;
    int64_t _childNamesOffset;
};

class MMapGenome : public Genome {
    friend class MMapDNAIterator;
public:
    MMapTopSegmentData *getTopSegmentPointer(hal_index_t index);
    MMapBottomSegmentData *getBottomSegmentPointer(hal_index_t index);

    const std::string& getName() const;

    void setDimensions(
        const std::vector<hal::Sequence::Info>& sequenceDimensions,
        bool storeDNAArrays);

    void updateTopDimensions(
        const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

    void updateBottomDimensions(
        const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

    hal_size_t getNumSequences() const;
   
    Sequence* getSequence(const std::string& name);

    const Sequence* getSequence(const std::string& name) const;

    Sequence* getSequenceBySite(hal_size_t position);
    const Sequence* getSequenceBySite(hal_size_t position) const;
   
    SequenceIteratorPtr getSequenceIterator(
        hal_index_t position);

    SequenceIteratorConstPtr getSequenceIterator(
        hal_index_t position) const;

    SequenceIteratorConstPtr getSequenceEndIterator() const;

    MetaData* getMetaData();

    const MetaData* getMetaData() const;

    Genome* getParent();

    const Genome* getParent() const;

    Genome* getChild(hal_size_t childIdx);

    const Genome* getChild(hal_size_t childIdx) const;

    hal_size_t getNumChildren() const;

    hal_index_t getChildIndex(const Genome* child) const;

    bool containsDNAArray() const;

    const Alignment* getAlignment() const;

    void rename(const std::string &newName);

    // SEGMENTED SEQUENCE INTERFACE

    hal_size_t getSequenceLength() const;
   
    hal_size_t getNumTopSegments() const;

    hal_size_t getNumBottomSegments() const;

    TopSegmentIteratorPtr getTopSegmentIterator(
        hal_index_t position);

    TopSegmentIteratorConstPtr getTopSegmentIterator(
        hal_index_t position) const;

    TopSegmentIteratorConstPtr getTopSegmentEndIterator() const;

    BottomSegmentIteratorPtr getBottomSegmentIterator(
        hal_index_t position);

    BottomSegmentIteratorConstPtr getBottomSegmentIterator(
        hal_index_t position) const;

    BottomSegmentIteratorConstPtr getBottomSegmentEndIterator() const;

    DNAIteratorPtr getDNAIterator(hal_index_t position);

    DNAIteratorConstPtr getDNAIterator(hal_index_t position) const;

    DNAIteratorConstPtr getDNAEndIterator() const;

    ColumnIteratorConstPtr getColumnIterator(const std::set<const Genome*>* targets,
                                             hal_size_t maxInsertLength,
                                             hal_index_t position,
                                             hal_index_t lastPosition,
                                             bool noDupes,
                                             bool noAncestors,
                                             bool reverseStrand,
                                             bool unique,
                                             bool onlyOrthologs) const;

    void getString(std::string& outString) const;

    void setString(const std::string& inString);

    void getSubString(std::string& outString, hal_size_t start,
                      hal_size_t length) const;

    void setSubString(const std::string& intString, 
                      hal_size_t start,
                      hal_size_t length);

    RearrangementPtr getRearrangement(hal_index_t position,
                                      hal_size_t gapLengthThreshold,
                                      double nThreshold,
                                      bool atomic = false) const;
   
    GappedTopSegmentIteratorConstPtr getGappedTopSegmentIterator(
        hal_index_t i, hal_size_t gapThreshold, bool atomic) const;

    GappedBottomSegmentIteratorConstPtr getGappedBottomSegmentIterator(
        hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
        bool atomic) const;

protected:
    char *getDNAArray() { return ((char *) _data) + _data->_dnaOffset; }

private:
    void setSequenceElement(size_t i, const std::string &name, hal_size_t startPos);
    MMapAlignment *_alignment;
    MMapGenomeData *_data;
    mutable std::map<hal_size_t, MMapSequence*> _sequencePosCache;
    mutable std::vector<MMapSequence*> _zeroLenPosCache;
    mutable std::map<std::string, MMapSequence*> _sequenceNameCache;
    void deleteSequenceCache();
    void loadSequencePosCache() const;
    void loadSequenceNameCache() const;
};
}
#endif
