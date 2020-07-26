/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halDefs.h"
#include "halSequence.h"
#include "halGenome.h"


/* thrown when sequence not found in genome */
hal::SequenceNotFoundException::SequenceNotFoundException(const Genome* genome,
                                                          const std::string &name):
    hal_exception("Sequence '" + name + "' not found in genome '" + genome->getName() + "'") {
}
