/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFSCANREFERENCE_H
#define _HALMAFSCANREFERENCE_H

#include "halMafScanner.h"
#include <cstdlib>
#include <deque>
#include <iostream>
#include <map>
#include <string>

namespace hal {

    /** Parse a MAF file line by line, getting some dimension stats
     * and maybe checking for some errros. */
    class MafScanReference : private MafScanner {
      public:
        MafScanReference();
        ~MafScanReference();

        std::string getRefName(const std::string &mafPath);

      private:
        void aLine();
        void sLine();
        void end();

        std::string _name;
    };
}

#endif
// Local Variables:
// mode: c++
// End:
