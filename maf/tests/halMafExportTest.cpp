/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafTests.h"
#include "halMafExport.h"

using namespace std;
using namespace hal;

CuSuite *halMafExportTestSuite(void)
{
  CuSuite* suite = CuSuiteNew();
  return suite;
}
