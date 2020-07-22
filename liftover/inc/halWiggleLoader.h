/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALWIGGLELOADER_H
#define _HALWIGGLELOADER_H

#include "halWiggleScanner.h"
#include "halWiggleTiles.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace hal {

    /** Quick hack to load a wiggle into memory, in order to have wiggleLiftover
     * --append option work better.  Ideally would be base class of WiggleLiftover
     * but don't have time to refactor right now */
    class WiggleLoader : public WiggleScanner {
      public:
        WiggleLoader();
        virtual ~WiggleLoader();

        void load(AlignmentConstPtr alignment, const Genome *genome, std::istream *inputFile, WiggleTiles<double> *vals);

      protected:
        virtual void visitLine();
        virtual void visitHeader();

        AlignmentConstPtr _alignment;
        const Genome *_srcGenome;
        const Sequence *_srcSequence;
        WiggleTiles<double> *_vals;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
