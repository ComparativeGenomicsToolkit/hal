/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTMAPPER_H
#define _HALSEGMENTMAPPER_H
#include "halDefs.h"
#include "halSegmentIterator.h"
#include <set>

namespace hal {
    class Segment;
    class MappedSegmentSet;
    class Genome;

    /** Get homologous segments in target genome.  Returns the number
      * of mapped segments found.
      * @param source Input.
      * @param outSegments  Output.  Mapped segments are sorted along the
      * *target* genome.
      * @param tgtGenome  Target genome to map to.  Can be the same as current.
      * @param genomesOnPath Intermediate genomes that must be visited
      * on the way down from coalescenceLimit to tgt.  If this is
      * specified as NULL, then the path will be computed automatically
      * (using hal::getGenomesInSpanningTree(coalescenceLimit, tgtGenome)).
      * Specifying this can avoid recomputing the path over and over again
      * when, say, calling halMapSegment repeatedly for the same
      * source and target.
      * @param doDupes  Specify whether paralogy edges are followed
      * @param minLength Minimum length of segments to consider.  It is
      * potentially much faster to filter using this parameter than
      * doing a second pass on the output.  If minLength is 0, then no
      * segments are filtered based on length.
      * @param coalescenceLimit Any paralogs that coalesce in or below
      * this genome will be mapped to the target as well. Must be the
      * MRCA or higher. By default, the coalescenceLimit is the MRCA.
      * @param mrca The MRCA of the source and target genomes. By
      * default, it is computed automatically. */
    hal_size_t halMapSegment(const SegmentIterator *source, MappedSegmentSet &outSegments, const Genome *tgtGenome,
                             const std::set<const Genome *> *genomesOnPath = NULL, bool doDupes = true,
                             hal_size_t minLength = 0, const Genome *coalescenceLimit = NULL, const Genome *mrca = NULL);

    /* call main function with smart pointer */
    hal_size_t halMapSegmentSP(const SegmentIteratorPtr &source, MappedSegmentSet &outSegments, const Genome *tgtGenome,
                               const std::set<const Genome *> *genomesOnPath = NULL, bool doDupes = true,
                               hal_size_t minLength = 0, const Genome *coalescenceLimit = NULL, const Genome *mrca = NULL);
}
#endif
