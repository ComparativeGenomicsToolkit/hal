/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _HDF5DNADRIVER_H
#define _HDF5DNADRIVER_H
#include "halDnaDriver.h"

namespace hal {
    class Hdf5Genome;
    class Hdf5ExternalArray;

    /**
     * HDF5 implementation of DnaAccess.
     */
    class HDF5DnaAccess : public DnaAccess {
      public:
        HDF5DnaAccess(Hdf5Genome *genome, Hdf5ExternalArray *dnaArray, hal_index_t index);

        virtual ~HDF5DnaAccess() {
        }

        void flush();

      protected:
        virtual void fetch(hal_index_t index) const;

      private:
        Hdf5ExternalArray *_dnaArray;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
