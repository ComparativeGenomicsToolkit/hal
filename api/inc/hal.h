/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ___HAL_H
#define ___HAL_H

/** \mainpage HAL API
 *
 * \section intro_sec Introduction
 *
 * This API is reads and writes genomic alignments in the 
 * Hierarchical ALignment (HAL) format.  The interface is independent
 * of how the alignment is stored.  All functionality is accessed
 * via an Alignment object as a starting point.  Alignment objects are
 * created with the XXXXAlignmentInstance() functions (where XXXX is
 * the implementation).  The HAL API is designed so that the user never
 * explicitly allocates or frees any memory associated with it (doing
 * so should in theory be prevented at compile-time). 
 *
 * \section hdf5_sec HDF5 Implementation
 *
 * Currently only a prototye (and still under developmenet)
 * HDF5 implentation is provided. It requires 
 * a current version of the HDF5 library 
 * (http://www.hdfgroup.org/downloads/index.html) to be installed.  
 *
 */

#include "counted_ptr.h"
#include "halDefs.h"
#include "halCommon.h"
#include "halAlignmentInstance.h"
#include "halAlignment.h"
#include "halGenome.h"
#include "halBottomSegment.h"
#include "halBottomSegmentIterator.h"
#include "halMetaData.h"
#include "halSequence.h"
#include "halSequenceIterator.h"
#include "halDNAIterator.h"
#include "halSegment.h"
#include "halSegmentIterator.h"
#include "halSegmentedSequence.h"
#include "halTopSegment.h"
#include "halTopSegmentIterator.h"

#endif
