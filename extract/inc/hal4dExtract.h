/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALEXTRACT4D_H
#define _HALEXTRACT4D_H

#include "halBedScanner.h"
#include <deque>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace hal {

    /** Use the halBedScanner to parse a bed file, extracting 4d codons from
     * each position */
    class Extract4d : public BedScanner {
      public:
        Extract4d();
        virtual ~Extract4d();

        void run(const Genome *refGenome, std::istream *inBedStream, std::ostream *outBedStream, bool conserved = false);

        static const char CodonPrefixTable[2][8];

      protected:
        virtual void visitLine();

        void extractBlocks4d(bool conserved);
        void write();

      protected:
        std::ostream *_outBedStream;
        const Genome *_refGenome;
        const Sequence *_refSequence;
        std::deque<BedLine> _outBedLines;
        bool _conserved;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
