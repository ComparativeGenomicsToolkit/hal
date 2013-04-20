/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBLOCKINTERPOLATE_H
#define _HALBLOCKINTERPOLATE_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"

namespace hal {

void blockInterpolateGenome(const Genome* genome, const Sequence* sequence,
                            const Genome* target,
                            hal_size_t start, hal_size_t length,
                            hal_size_t step);

void blockInterpolateSequence(const Sequence* sequence, const Genome* target,
                              hal_size_t start, hal_size_t length,
                              hal_size_t step);

}

#endif
