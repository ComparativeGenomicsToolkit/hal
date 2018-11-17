/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _HDF5STORAGEDRIVER_H
#define _HDF5STORAGEDRIVER_H
#include "halStorageDriver.h"

namespace hal {
    class Hdf5Genome;
    class Hdf5ExternalArray;
    
    /**
     * HDF5 implementation of DnaAccess.
     */
    class HDF5DnaAccess: public DnaAccess {
        public:
        HDF5DnaAccess(Hdf5Genome* genome,
                       Hdf5ExternalArray* dnaArray,
                       hal_index_t index);
        
        virtual ~HDF5DnaAccess() {
        }

        void flush();

        protected:
        virtual void fetch(hal_index_t index) const;

        private:
        Hdf5Genome* _genome;
        Hdf5ExternalArray* _dnaArray;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
