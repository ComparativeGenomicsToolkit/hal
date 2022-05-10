/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOMETADATA_H
#define _YOMOMETADATA_H

#include "halMetaData.h"
#include "yomoExternalArray.h"
#include <H5Cpp.h>
#include <map>
#include <string>

namespace hal {

    /**
     * YOMO string map used for general metadata
     */
    class YOMOMetaData : public MetaData {
      public:
        YOMOMetaData();
        YOMOMetaData(H5::PortableH5Location *parent, const std::string &name);
        virtual ~YOMOMetaData();

        void set(const std::string &key, const std::string &value);
        const std::string &get(const std::string &key) const;
        bool has(const std::string &key) const;
        const std::map<std::string, std::string> &getMap() const;

        void write();

        void open(H5::PortableH5Location *parent, const std::string &name);

      private:
        static const std::string MetaGroupName;

        H5::PortableH5Location *_parent;
        H5::Group _group;
        std::map<std::string, std::string> _map;
        bool _dirty;
        std::string _name;
    };
}
#endif

// Local Variables:
// mode: c++
// End:
