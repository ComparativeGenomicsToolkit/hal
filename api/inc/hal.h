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
 * \section coord_sec Note on Coordinates
 *
 * Internally, all coordinates are relative to the forward strand of the
 * genome (as opposed to sequence).  Currently, all interface methods
 * that return coordinates (ie SegmentIterator::getStartPosition) return
 * values in this coordinate system.  Iterators can be obtained relative
 * to sequence coordinates (using the sequence interface), and the
 * reverse compliment is taken into account when moving iterators vertically,
 * but coordinate values themselves are by default always expressed relative
 * to the entire genome in the forward strand. 
 *
 * \section pointer_sec Note on pointers
 *
 * When a genome is opened, there is a unique pointer to it that is valid
 * until it is explicitly closed (or the alignment is closed). This guarantee
 * does not hold for any other types (ie sequences or segments).
 * All pointers that are retrived via iterators have lifespans 
 * limited by their iterators.
 * Note that this means that the following code:
 * TopSegment* ts = genome->getTopSegmentIterator()->getTopSegment()
 * returns an invalid pointer, since the iterator gets descructed immediately.  
 * They are not
 * unique either, so multiple pointers / iterators can contain the same 
 * underlying object. 
 *
 * \section using_sec Using
 *
 * include hal/lib/hal.h and link to hal/lib/halLib.a
 *
 * \section hdf5_sec HDF5 Implementation
 *
 * Currently only a prototye (and still under developmenet)
 * HDF5 implentation is provided. It requires 
 * a current version of the HDF5 library 
 * (http://www.hdfgroup.org/downloads/index.html) to be installed.  
 * It must be built with h5c++ (as opposed to g++).  
 *
 */

#include "counted_ptr.h"
#include "halDefs.h"
#include "halCommon.h"
#include "halAlignmentInstance.h"
#include "halCLParserInstance.h"
#include "halAlignment.h"
#include "halCLParser.h"
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
#include "halDNAIterator.h"
#include "halValidate.h"
#include "halColumnIterator.h"

#endif
