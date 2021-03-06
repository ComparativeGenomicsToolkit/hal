/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMETADATA_H
#define _HALMETADATA_H

#include "halDefs.h"
#include <map>
#include <string>

namespace hal {

    /**
     * Interface for alignment (or genome) metadata
     * MetaData is a set of key/value pairs where each key and each
     * value is represented by a string.
     */
    class MetaData {
      public:
        /** Destructor */
        virtual ~MetaData() {
        }

        /** Get read-only reference to the map of metadata */
        virtual const std::map<std::string, std::string> &getMap() const = 0;

        /** Get the value associated with a key (throws error if key doesn't exist)
         * @param key MetaData key */
        virtual const std::string &get(const std::string &key) const = 0;

        /** Set a key-value pair (create's if doesn't exist, updates if does)
         * @param key Key to update
         * @param value Value to update */
        virtual void set(const std::string &key, const std::string &value) = 0;

        /** Determine if key exists in metadata
         * @param key Key to test */
        virtual bool has(const std::string &key) const = 0;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
