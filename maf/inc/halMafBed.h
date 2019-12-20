/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFBED_H
#define _HALMAFBED_H

#include "halBedScanner.h"
#include "halMafExport.h"
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace hal {

    /** Use the halBedScanner to parse a bed file, running mafExport on each
     * line */
    class MafBed : public BedScanner {
      public:
        MafBed(std::ostream &mafStream, AlignmentConstPtr alignment, const Genome *refGenome,
               std::set<const Genome *> &targetSet, MafExport &mafExport);
        virtual ~MafBed();

        void run(std::istream *bedStream);

      protected:
        virtual void visitLine();

      protected:
        std::ostream &_mafStream;
        AlignmentConstPtr _alignment;
        const Genome *_refGenome;
        const Sequence *_refSequence;
        hal_index_t _refStart;
        hal_size_t _refLength;
        std::set<const Genome *> &_targetSet;
        MafExport &_mafExport;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
