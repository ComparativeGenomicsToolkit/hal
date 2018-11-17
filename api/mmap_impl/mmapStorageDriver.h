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
     * Mmap implementation of DnaAccess.
     */
    class MMapDnaAccess: public DnaAccess {
        public:
        MMapDnaAccess(MMapGenome* genome,
                      hal_index_t index);
        
        virtual ~MMapDnaAccess() {
        }

        void flush();

        protected:
        virtual void fetch(hal_index_t index) const;

        private:
        MMapGenome* _genome;
        bool _isUdcProtocol;
    };

}

#endif
// Local Variables:
// mode: c++
// End:
