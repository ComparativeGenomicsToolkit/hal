/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMODNAARRAY_H
#define _YOMODNAARRAY_H

#include "rawH5ExternalArray.h"
#include <H5Cpp.h>
#include <cassert>

namespace hal {

    /**
     * Wraps the RawH5ExternalArray with interface tailored to storing and accessing
     * only DNA characters
     */
    class YomoDnaArray {
      public:
        /** Constructor */
        YomoDnaArray();

        /** Destructor */
        virtual ~YomoDnaArray();

        /** Create a new array (overloads method in parent)
         * @param file YOMO file in which to add new array dataset
         * @param path location of new array in file
         * @param size Fixed length of the new array
         * @param cparams Creation parameters for new array (chunking, zipping) */
        void create(H5File *file, const std::string &path, hsize_t size,
                    const H5::DSetCreatPropList &cparms = H5::DSetCreatPropList::DEFAULT);

        /** Open an existing array
         * @param file YOMO file containing array to open
         * @param path location of array in file */
        void open(H5File *file, const std::string &path);

        /** Write any unsaved buffer contents back to the file */
        void write();

        /** Get read/write iterator
         * @param offset position of iterator in array */
        YomoDnaIterator getDnaIterator(hsize_t offset = 0);

        /** Get read-only iterator
         * @param offset position of iterator in array */
        YomoDnaConstIterator getDnaConstIterator(hsize_t offset = 0);

        /** Get size of array */
        hsize_t size();

      private:
        RawH5ExternalArray _array;
    };

    // INLINE METHODS

    inline YomoDnaArray::getDnaIterator(hsize_t offset) {
        assert(offset < size());
        return DnaIterator(_array, offset);
    }

    inline YomoDnaConstIterator(hsize_t offset) {
        assert(offset < size());
        return DnaIterator(_array, offset);
    }
}
// Local Variables:
// mode: c++
// End:
