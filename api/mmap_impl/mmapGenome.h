#ifndef _MMAPGENOME_H
#define _MMAPGENOME_H
#include <map>
#include "halGenome.h"
#include "mmapAlignment.h"
#include "mmapString.h"
#include "mmapMetaData.h"
#include "mmapTopSegmentData.h"
#include "mmapBottomSegmentData.h"
#include "mmapPerfectHashTable.h"
#include "mmapGenomeSiteMap.h"

namespace hal {
class MMapBottomSegmentData;
class MMapSequence;
class MMapSequenceData;
class MmapAlignment;

class MMapGenomeData {
    friend class MMapGenome;
public:
    char *getDNA(MMapAlignment *alignment, size_t start, size_t length) const;
    std::string getName(MMapAlignment *alignment) const;
    void setName(MMapAlignment *alignment, const std::string &name);
    void initializeName(MMapAlignment *alignment, const std::string &name);
    MMapTopSegmentData *getTopSegmentData(MMapAlignment *alignment, hal_index_t index);
    MMapBottomSegmentData *getBottomSegmentData(MMapAlignment *alignment, MMapGenome *genome, hal_index_t index);
private:
    hal_size_t _totalSequenceLength;
    hal_size_t _numSequences;
    hal_size_t _numMetadata;
    hal_size_t _numTopSegments;
    hal_size_t _numBottomSegments;

    size_t _nameOffset;
    size_t _sequenceHashOffset;
    size_t _genomeSiteMapOffset;
    size_t _sequencesOffset;
    size_t _metadataOffset;
    size_t _dnaOffset;
    size_t _topSegmentsOffset;
    size_t _bottomSegmentsOffset;
    size_t _childNamesOffset;
};

class MMapGenome : public Genome {
public:
    MMapGenome(MMapAlignment *alignment, MMapGenomeData *data, size_t arrayIndex) :
        Genome(alignment, data->getName(alignment)),
        _alignment(alignment),
        _data(data),
        _arrayIndex(arrayIndex),
        _name(data->getName(_alignment)),
        _metaData(_alignment, _data->_metadataOffset),
        _sequenceNameHash(alignment->getMMapFile(), data->_sequenceHashOffset),
        _genomeSiteMap(alignment->getMMapFile(), data->_genomeSiteMapOffset) {
        _sequenceObjCache.resize(data->_numSequences);
    };
    MMapGenome(MMapAlignment *alignment, MMapGenomeData *data, size_t arrayIndex, const std::string &name) :
        Genome(alignment, name),
        _alignment(alignment),
        _data(data),
        _arrayIndex(arrayIndex),
        _name(name),
        _metaData(_alignment),
        _sequenceNameHash(alignment->getMMapFile(), data->_sequenceHashOffset),
        _genomeSiteMap(alignment->getMMapFile(), data->_genomeSiteMapOffset) {
        _data->initializeName(_alignment, _name);
        _data->_metadataOffset = _metaData.getOffset();
        _sequenceObjCache.resize(data->_numSequences);
    };

    virtual ~MMapGenome();

    MMapTopSegmentData *getTopSegmentPointer(hal_index_t index) { return _data->getTopSegmentData(_alignment, index); };
    MMapBottomSegmentData *getBottomSegmentPointer(hal_index_t index) { return _data->getBottomSegmentData(_alignment, this, index);  };

    void updateGenomeArrayBasePtr(MMapGenomeData *base) {
        _data = base + _arrayIndex;
    }

    const std::string& getName() const;

    hal_index_t getArrayIndex() const {
        return _arrayIndex;
    }

    void setDimensions(
        const std::vector<hal::Sequence::Info>& sequenceDimensions,
        bool storeDNAArrays);

    void updateTopDimensions(
        const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

    void updateBottomDimensions(
        const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions);

    hal_size_t getNumSequences() const;

    // FIXME: these should be in interface and hdf5
    Sequence* getSequenceByIndex(hal_index_t index);

    const Sequence* getSequenceByIndex(hal_index_t index) const;
    
    Sequence* getSequence(const std::string& name);

    const Sequence* getSequence(const std::string& name) const;

    Sequence* getSequenceBySite(hal_size_t position);
    const Sequence* getSequenceBySite(hal_size_t position) const;
   
    SequenceIteratorPtr getSequenceIterator(
        hal_index_t position);

    SequenceIteratorPtr getSequenceIterator(
        hal_index_t position) const;

    MetaData* getMetaData();

    const MetaData* getMetaData() const;

    bool containsDNAArray() const;

    const Alignment* getAlignment() const;   // can't be inlined due to mutual include

    Alignment* getAlignment();  // can't be inlined due to mutual include

    void rename(const std::string &newName);

    // SEGMENTED SEQUENCE INTERFACE

    hal_size_t getSequenceLength() const;
   
    hal_size_t getNumTopSegments() const;

    hal_size_t getNumBottomSegments() const;

    TopSegmentIteratorPtr getTopSegmentIterator(
        hal_index_t position);

    TopSegmentIteratorPtr getTopSegmentIterator(
        hal_index_t position) const;

    BottomSegmentIteratorPtr getBottomSegmentIterator(
        hal_index_t position);

    BottomSegmentIteratorPtr getBottomSegmentIterator(
        hal_index_t position) const;

    DnaIteratorPtr getDnaIterator(hal_index_t position);

    DnaIteratorPtr getDnaIterator(hal_index_t position) const;

    ColumnIteratorPtr getColumnIterator(const std::set<const Genome*>* targets,
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
   
    GappedTopSegmentIteratorPtr getGappedTopSegmentIterator(
        hal_index_t i, hal_size_t gapThreshold, bool atomic) const;

    GappedBottomSegmentIteratorPtr getGappedBottomSegmentIterator(
        hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
        bool atomic) const;

    MMapAlignment *_alignment;
    MMapSequenceData *getSequenceData(size_t i) const;

    char *getDNA(size_t start, size_t length) {
        return _data->getDNA(_alignment, start, length);
    }

private:
    void createSequenceNameHash(size_t numSequences);
    void createGenomeSiteMap(size_t numSequences);
    void setSequenceData(size_t i, hal_index_t startPos, hal_index_t topSegmentStartIndex,
                         hal_index_t bottomSegmentStartIndex, const Sequence::Info &sequenceInfo);
    std::vector<Sequence::UpdateInfo> getCompleteInputDimensions(const std::vector<Sequence::UpdateInfo>& inputDimensions, bool isTop);
    void deleteSequenceCache();

    MMapGenomeData *_data;
    size_t _arrayIndex;     // Index within the alignment's genome array.
    std::string _name;
    MMapMetaData _metaData;
    MMapPerfectHashTable _sequenceNameHash;
    MMapGenomeSiteMap _genomeSiteMap;
    
    mutable std::vector<MMapSequence*> _sequenceObjCache;
};

inline std::string MMapGenomeData::getName(MMapAlignment *alignment) const {
    return MMapString(alignment, _nameOffset).c_str();
}

inline void MMapGenomeData::initializeName(MMapAlignment *alignment, const std::string &nameStr) {
    MMapString name{alignment, nameStr};
    _nameOffset = name.getOffset();
}

inline void MMapGenomeData::setName(MMapAlignment *alignment, const std::string &newName) {
    MMapString name{alignment, _nameOffset};
    _nameOffset = name.set(newName);
}

inline MMapTopSegmentData *MMapGenomeData::getTopSegmentData(MMapAlignment *alignment, hal_index_t index) {
    // We request twice the segment length here because checking the length of
    // this segment requires reading the start position of the following
    // segment.
    return static_cast<MMapTopSegmentData*>(alignment->resolveOffset(_topSegmentsOffset + index * sizeof(MMapTopSegmentData), 2 *sizeof(MMapTopSegmentData)));
}

inline MMapBottomSegmentData *MMapGenomeData::getBottomSegmentData(MMapAlignment *alignment, MMapGenome *genome, hal_index_t index) {
    size_t segmentSize = MMapBottomSegmentData::getSize(genome);
    // We request twice the segment length here because checking the length of
    // this segment requires reading the start position of the following
    // segment.
    return static_cast<MMapBottomSegmentData*>(alignment->resolveOffset(_bottomSegmentsOffset + index * segmentSize, 2 * segmentSize));
}

inline char *MMapGenomeData::getDNA(MMapAlignment *alignment, size_t start, size_t length) const {
    return static_cast<char*>(alignment->resolveOffset(_dnaOffset + start, length));
}
}
#endif
// Local Variables:
// mode: c++
// End:
