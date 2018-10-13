#ifndef _MMAPGENOME_H
#define _MMAPGENOME_H
#include <map>
#include "halGenome.h"
#include "mmapAlignment.h"
namespace hal {
class MMapTopSegmentData;
class MMapBottomSegmentData;
class MMapSequence;
class MMapSequenceData;

class MMapGenomeData {
    friend class MMapGenome;
public:
    std::string getName(MMapAlignment *alignment) const;
    void setName(MMapAlignment *alignment, const std::string &name);
protected:
    hal_size_t _nameLength;
    hal_size_t _totalSequenceLength;
    hal_size_t _numSequences;
    hal_size_t _numMetadata;
    hal_size_t _numTopSegments;
    hal_size_t _numBottomSegments;
    hal_size_t _numChildren;

    size_t _nameOffset;
    size_t _sequencesOffset;
    size_t _metadataOffset;
    size_t _dnaOffset;
    size_t _topSegmentsOffset;
    size_t _bottomSegmentsOffset;
    size_t _childNamesOffset;
};

class MMapGenome : public Genome {
    friend class MMapDNAIterator;
public:
    MMapGenome(MMapAlignment *alignment, MMapGenomeData *data) : _alignment(alignment), _data(data) {
        _name = _data->getName(_alignment);
    };
    MMapGenome(MMapAlignment *alignment, MMapGenomeData *data, const std::string &name) :
        _alignment(alignment), _data(data), _name(name) {
        _data->setName(_alignment, _name);
    };

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

    MMapAlignment *_alignment;
    MMapSequenceData *getSequenceData(size_t i) const;
protected:
    char *getDNAArray() { return ((char *) _data) + _data->_dnaOffset; }

private:
    void setSequenceData(size_t i, hal_index_t startPos, hal_index_t topSegmentStartIndex,
                         hal_index_t bottomSegmentStartIndex, const Sequence::Info &sequenceInfo);
    MMapGenomeData *_data;
    std::string _name;
    mutable std::map<hal_size_t, MMapSequence*> _sequencePosCache;
    mutable std::vector<MMapSequence*> _zeroLenPosCache;
    mutable std::map<std::string, MMapSequence*> _sequenceNameCache;
    void deleteSequenceCache();
    void loadSequencePosCache() const;
    void loadSequenceNameCache() const;
};

std::string MMapGenomeData::getName(MMapAlignment *alignment) const {
    return (const char *) alignment->resolveOffset(_nameOffset, _nameLength);
};
void MMapGenomeData::setName(MMapAlignment *alignment, const std::string &newName) {
    size_t size = newName.size() + 1;
    _nameOffset = alignment->allocateNewArray(sizeof(char) * size);
    strncpy((char *) alignment->resolveOffset(_nameOffset, size), newName.c_str(), size);
    _nameLength = size;
};
}
#endif
