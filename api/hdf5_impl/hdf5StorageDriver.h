/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _HDF5STORAGEDRIVER_H
#define _HDF5STORAGEDRIVER_H
#include "halStorageDriver.h"

namespace hal {
    class HDF5Genome;
    class HDF5ExternalArray;
    
    /**
     * HDF5 implementation of DNAAccess.
     */
    class HDF5DNAAccess: public DNAAccess {
        public:
        HDF5DNAAccess(HDF5Genome* genome,
                       HDF5ExternalArray* dnaArray,
                       hal_index_t index);
        
        virtual ~HDF5DNAAccess() {
        }

        void flush();

        protected:
        virtual void fetch(hal_index_t index) const;

        private:
        HDF5Genome* _genome;
        HDF5ExternalArray* _dnaArray;
    };
}

#endif
