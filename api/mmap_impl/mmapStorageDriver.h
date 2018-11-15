/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _MMAPSTORAGEDRIVER_H
#define _MMAPSTORAGEDRIVER_H
#include "halStorageDriver.h"

namespace hal {
    class MMapGenome;
    class MMapAlignment;
    
    /**
     * Mmap implementation of DNAStorage.
     */
    class MMapDNAStorage: public DNAStorage {
        public:
        MMapDNAStorage(MMapGenome* genome,
                       hal_index_t index);
        
        virtual ~MMapDNAStorage() {
        }

        protected:
        virtual void fetch(hal_index_t index) const;

        private:
        MMapGenome* _genome;
        bool _isUdcProtocol;
    };

}

#endif
