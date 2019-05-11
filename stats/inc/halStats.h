/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSTATS_H
#define _HALSTATS_H

#include "hal.h"
#include <iostream>
#include <string>
#include <vector>

namespace hal {

    struct GenomeStats : public hal::Sequence::Info {
        size_t _numChildren;
        size_t _numSequences;
    };

    class HalStats {
      public:
        HalStats();
        HalStats(const Alignment *alignment);
        virtual ~HalStats();

        void printCsv(std::ostream &outStream) const;
        void readAlignment(const Alignment *alignment);

      protected:
        void readGenomeRecursive(const Alignment *alignment, const Genome *genome);

        std::string _tree;
        std::vector<GenomeStats> _genomeStatsVec;
    };
}

inline std::ostream &operator<<(std::ostream &os, const hal::HalStats &halStats) {
    halStats.printCsv(os);
    return os;
}

#endif
// Local Variables:
// mode: c++
// End:
