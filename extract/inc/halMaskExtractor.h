/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMASKEXTRACTOR_H
#define _HALMASKEXTRACTOR_H

#include "hal.h"
#include <iostream>
#include <string>
#include <vector>

namespace hal {

    class MaskExtractor {
      public:
        MaskExtractor();
        virtual ~MaskExtractor();

        void extract(const Alignment *alignment, const Genome *genome, std::ostream *bedStream, hal_size_t extend,
                     double extendPct);

      protected:
        void addMaskedBasesToCache();
        void extendCachedIntervals();
        void writeCachedIntervals();

      protected:
        AlignmentConstPtr _alignment;
        const Genome *_genome;
        const Sequence *_sequence;
        std::ostream *_bedStream;
        hal_size_t _extend;
        double _extendPct;
        PositionCache _posCache;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
