#ifndef _MMAPSEQUENCE_H
#define _MMAPSEQUENCE_H
#include "halSequence.h"
namespace hal {
class MMapSequence : public Sequence {
public:
   MMapSequence(MMapGenome* genome,
                hal_index_t index) {};

   // SEQUENCE INTERFACE
   std::string getName() const;

   std::string getFullName() const;

   const Genome* getGenome() const;

   Genome* getGenome();

   hal_index_t getStartPosition() const;

   hal_index_t getEndPosition() const;

   hal_index_t getArrayIndex() const;

   hal_index_t getTopSegmentArrayIndex() const;

   hal_index_t getBottomSegmentArrayIndex() const;

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

   void setName(const std::string &newName);
};
}
#endif
