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
 * created with the openAlignment() and AlignmentInstance() functions 
 * found in api/inc/halAlignmentInstance.h
 * The HAL API is designed so that the user never
 * explicitly allocates or frees any memory associated with it. New 
 * iterators and such are created and returned inside reference-counted
 * pointers.  Other, more persistent objects such as genomes and sequences
 * are accessed via const pointers and are created and cleaned up through the
 * alignment class.
 *
 * Several tools are included in the API and can serve as examples for all
 * the difference classes.
 *
 * Please consult hal/README.md (or pdf) for more information on the HAL
 * format and tools. 
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
 * to the entire genome in the forward strand.  In other words, reversing an
 * iterator will reverse the start and end coordinates, but these coordinates
 * will remain relative to the first base of the genome in the forward strand. 
 *
 * \section pointer_sec Note on pointers
 *
 * When a genome is opened, there is a unique pointer to it that is valid
 * until it is explicitly closed (or the alignment is closed). This guarantee
 * does not hold for any other types (ie sequences or segments).
 * All pointers that are retrived via iterators have lifespans 
 * limited by their iterators (retrieving segment pointers from the iterators
 * is now deprecated and shouldn't be done so this is no longer an issue).
 *
 * \section using_sec Using
 *
 * include hal/lib/hal.h and link to hal/lib/halLib.a
 *
 */

#include "halDefs.h"
#include "halCommon.h"
#include "halPositionCache.h"
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
#include "halSlicedSegment.h"
#include "halMappedSegment.h"
#include "halSegmentIterator.h"
#include "halSegmentedSequence.h"
#include "halTopSegment.h"
#include "halTopSegmentIterator.h"
#include "halDNAIterator.h"
#include "halValidate.h"
#include "halColumnIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halRearrangement.h"

#endif
