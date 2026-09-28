/*
 * HAL format 3 API: everything a program reading an alignment needs.
 *
 *   auto aln = hal3::api::Alignment::open("apes.hal");
 *   hal3::api::Lifter lifter(aln);               // one per thread
 *   lifter.lift(aln->genome("hg38"), seq, start, end, aln->genome("Anc0"), true, pieces);
 *
 * Link libHal3.a, liftmap and twobit64 (and zlib, pthread; libcurl if liftmap was built
 * with HTTP).  C++11.  The C++ interface may change between releases; there is no ABI
 * promise.
 */
#ifndef _HAL3_API_H
#define _HAL3_API_H

#include "hal3Error.h"
#include "hal3Lift.h"
#include "hal3Reader.h"
#include "hal3Walk.h"

#endif
