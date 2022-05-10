/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOSEQUENCE_H
#define _YOMOSEQUENCE_H

#include "halSequence.h"
#include "yomoExternalArray.h"
#include "yomoGenome.h"
#include <H5Cpp.h>

namespace hal {

    class YomoSequenceIterator;

    class YomoSequence : public Sequence {
        friend class YomoSequenceIterator;

      public:
        YomoSequence(YomoGenome *genome, YomoExternalArray *idxArray, YomoExternalArray *nameArray, hal_index_t index);

        /** Destructor */
        ~YomoSequence();

        // SEQUENCE INTERFACE
        std::string getName() const;

        std::string getFullName() const;

        const Genome *getGenome() const;

        Genome *getGenome();

        hal_index_t getStartPosition() const;

        hal_index_t getEndPosition() const;

        hal_index_t getArrayIndex() const;

        hal_index_t getTopSegmentArrayIndex() const;

        hal_index_t getBottomSegmentArrayIndex() const;

        // SEGMENTED SEQUENCE INTERFACE

        hal_size_t getSequenceLength() const;

        hal_size_t getNumTopSegments() const;

        hal_size_t getNumBottomSegments() const;

        TopSegmentIteratorPtr getTopSegmentIterator(hal_index_t position);

        TopSegmentIteratorPtr getTopSegmentIterator(hal_index_t position) const;

        BottomSegmentIteratorPtr getBottomSegmentIterator(hal_index_t position);

        BottomSegmentIteratorPtr getBottomSegmentIterator(hal_index_t position) const;

        DnaIteratorPtr getDnaIterator(hal_index_t position);

        DnaIteratorPtr getDnaIterator(hal_index_t position) const;

        ColumnIteratorPtr getColumnIterator(const std::set<const Genome *> *targets, hal_size_t maxInsertLength,
                                            hal_index_t position, hal_index_t lastPosition, bool noDupes, bool noAncestors,
                                            bool reverseStrand, bool unique, bool onlyOrthologs) const;

        void getString(std::string &outString) const;

        void setString(const std::string &inString);

        void getSubString(std::string &outString, hal_size_t start, hal_size_t length) const;

        void setSubString(const std::string &intString, hal_size_t start, hal_size_t length);

        RearrangementPtr getRearrangement(hal_index_t position, hal_size_t gapLengthThreshold, double nThreshold,
                                          bool atomic = false) const;

        GappedTopSegmentIteratorPtr getGappedTopSegmentIterator(hal_index_t i, hal_size_t gapThreshold, bool atomic) const;

        GappedBottomSegmentIteratorPtr getGappedBottomSegmentIterator(hal_index_t i, hal_size_t childIdx,
                                                                      hal_size_t gapThreshold, bool atomic) const;

        void setName(const std::string &newName);

        // LOCAL NON-INTERFACE METHODS

        static H5::CompType idxDataType();
        static H5::StrType nameDataType(hal_size_t maxNameLength);

        void set(hal_size_t startPosition, const Sequence::Info &sequenceInfo, hal_size_t topSegmentStartIndex,
                 hal_size_t bottomSegmentStartIndex);

        void setNumTopSegments(hal_size_t numTopSegments);

        void setNumBottomSegments(hal_size_t numBottomSegments);

        void setTopSegmentArrayIndex(hal_size_t topIndex);

        void setBottomSegmentArrayIndex(hal_size_t bottomIndex);

      private:
        void refreshNameCache() const;
        void refreshFullNameCache() const;
        void refreshIdxCache() const;

        static const size_t startOffset;
        static const size_t topSegmentArrayIndexOffset;
        static const size_t bottomSegmentArrayIndexOffset;
        static const size_t totalSize;

        YomoExternalArray *_idxArray;
        YomoExternalArray *_nameArray;
        hal_index_t _index;
        YomoGenome *_genome;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
