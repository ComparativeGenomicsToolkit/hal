/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFEXPORT_H
#define _HALMAFEXPORT_H

#include "halMafBlock.h"
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace hal {

    class MafExport {
      public:
        MafExport():
            _mafStream(NULL), _maxRefGap(0), _noDupes(false), _noAncestors(false),
            _ucscNames(false), _unique(false), _append(false), _printTree(false),
            _onlyOrthologs(false), _keepEmptyRefBlocks(false) {
        }
        
        virtual ~MafExport() {
        }

        void convertSequence(std::ostream &mafStream, AlignmentConstPtr alignment, const Sequence *seq,
                             hal_index_t startPosition, hal_size_t length, const std::set<const Genome *> &targets);

        // Convert all columns in the leaf genomes to MAF. Each column is
        // reported exactly once regardless of the unique setting, although
        // this may change in the future. Likewise, maxRefGap has no
        // effect, although noDupes will work.
        void convertEntireAlignment(std::ostream &mafStream, AlignmentConstPtr alignment);

        void setMaxRefGap(hal_size_t maxRefGap) {
            _maxRefGap = maxRefGap;
        }
        void setNoDupes(bool noDupes) {
            _noDupes = noDupes;
        }
        void setNoAncestors(bool noAncestors) {
            _noAncestors = noAncestors;
        }
        void setUcscNames(bool ucscNames) {
            _ucscNames = ucscNames;
        }
        void setUnique(bool unique) {
            _unique = unique;
        }
        void setAppend(bool append) {
            _append = append;
        }
        void setMaxBlockLength(hal_index_t maxLength) {
            _mafBlock.setMaxLength(maxLength);
        }
        void setPrintTree(bool printTree) {
            _printTree = printTree;
        }
        void setOnlyOrthologs(bool onlyOrthologs) {
            _onlyOrthologs = onlyOrthologs;
        }
        void setKeepEmptyRefBlocks(bool keepEmptyRefBlocks) {
            _keepEmptyRefBlocks = keepEmptyRefBlocks;
        }

      protected:
        void writeHeader();

      protected:
        AlignmentConstPtr _alignment;
        std::ostream *_mafStream;
        MafBlock _mafBlock;
        hal_size_t _maxRefGap;
        bool _noDupes;
        bool _noAncestors;
        bool _ucscNames;
        bool _unique;
        bool _append;
        bool _printTree;
        bool _onlyOrthologs;
        bool _keepEmptyRefBlocks;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
