/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBLOCKLIFTOVER_H
#define _HALBLOCKLIFTOVER_H

#include "halLiftover.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace hal {

    class BlockLiftover : public Liftover {
      public:
        BlockLiftover();
        virtual ~BlockLiftover();

        /** Enable closing intermediate genomes during batch liftover
         *  to reduce peak memory.  Only safe when the caller does not
         *  hold open Genome pointers to intermediate genomes. */
        void setCloseGenomes(bool close) {
            _closeGenomes = close;
        }

      protected:
        void liftInterval(BedList &mappedBedLines);
        void liftIntervalBatchMap(hal_index_t globalStart, hal_index_t globalEnd, bool flip);
        void visitBegin();

        void cleanTargetParalogies();
        void readPSLInfo(std::vector<MappedSegmentPtr> &fragments, BedLine &outBedLine);

      protected:
        MappedSegmentSet _mappedSegments;
        SegmentIteratorPtr _refSeg;
        hal_index_t _lastIndex;
        std::set<const Genome *> _downwardPath;
        const Genome *_mrca;
        bool _useBatchPath;
        bool _closeGenomes;
        std::string _mrcaName;
        std::set<std::string> _downwardPathNames;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
