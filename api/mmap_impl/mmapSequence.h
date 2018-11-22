#ifndef _MMAPSEQUENCE_H
#define _MMAPSEQUENCE_H
#include "halSequence.h"
#include "mmapGenome.h"
#include "mmapSequenceData.h"

namespace hal {
class MMapSequence : public Sequence {
    friend class MMapSequenceIterator;
public:
    MMapSequence(MMapGenome *genome,
                 MMapSequenceData *data) : _genome(genome), _data(data) {}

    MMapSequence(MMapGenome *genome,
                 MMapSequenceData *data,
                 hal_index_t index,
                 hal_index_t startPosition,
                 hal_size_t length,
                 hal_index_t topSegmentStartIndex,
                 hal_index_t bottomSegmentStartIndex,
                 hal_size_t numTopSegments,
                 hal_size_t numBottomSegments,
                 const std::string &name) :
        _genome(genome),
        _data(data) {
        _data->_index = index;
        _data->_startPosition = startPosition;
        _data->_length = length;
        _data->_topSegmentStartIndex = topSegmentStartIndex;
        _data->_bottomSegmentStartIndex = bottomSegmentStartIndex;
        _data->_numTopSegments = numTopSegments;
        _data->_numBottomSegments = numBottomSegments;
        _data->setName(_genome->_alignment, name);
    };

   // SEQUENCE INTERFACE
   std::string getName() const { return _data->getName(_genome->_alignment); }

   std::string getFullName() const { return _genome->getName() + "." + _data->getName(_genome->_alignment); }

   const Genome* getGenome() const { return _genome; }

   Genome* getGenome() { return _genome; }

   hal_index_t getStartPosition() const { return _data->_startPosition; }

   hal_index_t getEndPosition() const { return _data->_startPosition + _data->_length; }

   hal_index_t getArrayIndex() const { return _data->_index; };

   hal_index_t getTopSegmentArrayIndex() const { return _data->_topSegmentStartIndex; };

   hal_index_t getBottomSegmentArrayIndex() const { return _data->_bottomSegmentStartIndex; };

   // SEGMENTED SEQUENCE INTERFACE

   hal_size_t getSequenceLength() const { return _data->_length; }
   
   hal_size_t getNumTopSegments() const { return _data->_numTopSegments; }

   hal_size_t getNumBottomSegments() const { return _data->_numBottomSegments; }

   void setNumTopSegments(hal_size_t topSegments) { _data->_numTopSegments = topSegments; }

   void setNumBottomSegments(hal_size_t bottomSegments) { _data->_numBottomSegments = bottomSegments; }

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

   void setName(const std::string &newName) { _data->setName(_genome->_alignment, newName); }

    void setTopSegmentStartIndex(hal_index_t index) { _data->_topSegmentStartIndex = index; }
    void setBottomSegmentStartIndex(hal_index_t index) { _data->_bottomSegmentStartIndex = index; }

    private:
    MMapGenome *_genome;
    MMapSequenceData *_data;
};
}
#endif
// Local Variables:
// mode: c++
// End:
