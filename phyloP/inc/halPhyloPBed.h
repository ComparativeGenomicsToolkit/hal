/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALPHYLOPBED_H
#define _HALPHYLOPBED_H

#include "halBedScanner.h"
#include "halPhyloP.h"
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace hal {

    /** Use the halBedScanner to parse a bed file, running halPhyloP on each
     * line */
    class PhyloPBed : public BedScanner {
      public:
        PhyloPBed(AlignmentConstPtr alignment, const Genome *refGenome, const Sequence *refSequence, hal_index_t start,
                  hal_size_t length, hal_size_t step, PhyloP &phyloP, std::ostream &outStream);
        virtual ~PhyloPBed();

        void run(std::istream *bedStream, int bedVersion = -1);

      protected:
        virtual void visitLine();

      protected:
        AlignmentConstPtr _alignment;
        const Genome *_refGenome;
        const Sequence *_refSequence;
        hal_index_t _refStart;
        hal_index_t _refLength;
        hal_size_t _step;
        PhyloP &_phyloP;
        std::ostream &_outStream;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
