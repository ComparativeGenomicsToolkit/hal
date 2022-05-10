/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _YOMOCOMMON_H
#define _YOMOCOMMON_H
#include <H5Cpp.h>

/* Class to disable YOMO exception printing.  It is re-enabled when the
 * class is destroyed.  Can comment out dontPrint() call for debugging.
 */
class YOMODisableExceptionPrinting {
  public:
    YOMODisableExceptionPrinting() {
        H5::Exception::getAutoPrint(_func, &_clientData);
        H5::Exception::dontPrint();
    }
    ~YOMODisableExceptionPrinting() {
        H5::Exception::setAutoPrint(_func, _clientData);
    }

  private:
    H5E_auto2_t _func;
    void *_clientData;
};

#endif

// Local Variables:
// mode: c++
// End:
