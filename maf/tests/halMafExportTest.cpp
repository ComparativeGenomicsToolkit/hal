/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafExport.h"
#include "halMafTests.h"

using namespace std;
using namespace hal;

CuSuite *halMafExportTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    return suite;
}
