/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALVALIDATE_H
#define _HALVALIDATE_H

#include <map>
#include <string>
#include <vector>
#include "halAlignment.h"

namespace hal {

/** Go through a bottom segment, and throw an exception if anything 
 * appears out of whack. */
void validateBottomSegment(const BottomSegment* bottomSegment);

/** Go through a top segment, and throw an exception if anything 
 * appears out of whack. */
void validateTopSegment(const TopSegment* topSegment);

/** Go through a sequence, and throw an exception if anything 
 * appears out of whack. */
void validateSequence(const Sequence* sequence);

/** Go through a genome, and throw an exception if anything 
 * appears out of whack. */
void validateGenome(const Genome* genome);

/** Go through a genome, and throw an exception if any duplications 
 * appears out of whack. */
void validateDuplications(const Genome* genome);

/** Go through an alignment, and throw an excpetion if anything
 * appears out of whack. */
void validateAlignment(AlignmentConstPtr alignment);

}
#endif

