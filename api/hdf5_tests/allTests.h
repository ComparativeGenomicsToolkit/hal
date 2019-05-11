/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _ALLTESTS_H
#define _ALLTESTS_H

extern "C" {
#include "CuTest.h"
}

CuSuite *hdf5TestSuite();
CuSuite *hdf5ExternalArrayTestSuite();
CuSuite *hdf5DNATypeTestSuite();
CuSuite *hdf5SegmentTypeTestSuite();
CuSuite *hdf5SequenceTypeTestSuite();

#endif
// Local Variables:
// mode: c++
// End:
