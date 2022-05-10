/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _YOMODNADRIVER_H
#define _YOMODNADRIVER_H
#include "halDnaDriver.h"

namespace hal {
    class YomoGenome;
    class YomoExternalArray;

    /**
     * YOMO implementation of DnaAccess.
     */
    class YOMODnaAccess : public DnaAccess {
      public:
        YOMODnaAccess(YomoGenome *genome, YomoExternalArray *dnaArray, hal_index_t index);

        virtual ~YOMODnaAccess() {
        }

        void flush();

      protected:
        virtual void fetch(hal_index_t index) const;

      private:
        YomoExternalArray *_dnaArray;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
