/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENT_H
#define _HALSEGMENT_H

#include <string>
#include <vector>
#include <set>
#include "halDefs.h"

namespace hal {

/** 
 * Interface for a segment of DNA. Note that segments should
 * not be written to outside of creating new genomes.
 */
class Segment
{
public:

   /** Set the current array index of the segment.  This writes no information
    * to the database, but just moves the position of the segment
    * @param genome Genome whose array we want to move segment to
    * @param arrayIndex Index in genomes array (in genome segment coordinates) 
    */
   virtual void setArrayIndex(Genome* genome, hal_index_t arrayIndex) = 0;

   /** Set the current array index of the segment.  This writes no information
    * to the database, but just moves the position of the segment
    * @param genome Genome whose array we want to move segment to
    * @param arrayIndex Index in genomes array (in genome segment coordinates) 
    */
   virtual void setArrayIndex(const Genome* genome, hal_index_t arrayIndex) 
   const = 0;
 
   /** Get the containing (read-only) genome */
   virtual const Genome* getGenome() const = 0;

   /** Get the containing genome */
   virtual Genome* getGenome() = 0;

   /** Get the containing (read-only) sequence */
   virtual const Sequence* getSequence() const = 0;

   /** Get the containing sequence */
   virtual Sequence* getSequence() = 0;

   /** Get the start position of the segment. Note that the start position
    * is currently always in FORWARD GENOME coordinates.  That said,
    * if the underlying object is an iterator is in reverse orientation, 
    * the segment is read from right (startposition) to left 
    * (startposition - length -1)
    * Similarly, if called from an iterator, the slicing offsets are taken
    * into account */
   virtual hal_index_t getStartPosition() const = 0;

   /** Get the segment's end position (start + len - 1) in the genome
    * Slicing and reversal will apply as above for iterators */
   virtual hal_index_t getEndPosition() const = 0;

   /** Get the length of the segment (number of bases)*/
   virtual hal_size_t getLength() const = 0;
   
   /** Get the DNA string corresponding to the segment from the genome 
    * @param outString string into which the results are copied */
   virtual void getString(std::string& outString) const = 0;

   /** Set the segment's start position in the genome 
    * @param startPos Start position 
    * @param length Length*/
   virtual void setCoordinates(hal_index_t startPos, hal_size_t length) = 0;

   /** Get the index of the segment in the segment array */
   virtual hal_index_t getArrayIndex() const = 0;

   /** Determine if current segment is to the left of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool leftOf(hal_index_t genomePos) const = 0;

   /** Determine if current segment is to the right of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool rightOf(hal_index_t genomePos) const = 0;

   /** Determine if current segment is to the right of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool overlaps(hal_index_t genomePos) const = 0;

   /** Check whether segment is the first segment of a sequence 
    * if underlying object is reversed iterator, then checks if last*/
   virtual bool isFirst() const = 0;

   /** Check whether segment is the last segment of a sequence
    * if underlying object is reversed iterator, then checks if first*/
   virtual bool isLast() const = 0;
   
   /** Test if the fraction of N's in the segment is greater than 
    * a given threshold.
    * @param nThreshold Maximum fraction of N's in the segment for it 
    * to no be considered missing data */
   virtual bool isMissingData(double nThreshold) const = 0;

   /** Check if underlying segment is a top segment (easier than doing a
    * downcast.  Returns true if it's a top segment and false if it's a 
    * bottom segment */
   virtual bool isTop() const = 0;

   /** Get homologous segments in target genome.  Returns the number
    * of mapped segments found.
    * @param outSegments  Output.  Mapped segments are sorted along the
    * *target* genome.
    * @param tgtGenome  Target genome to map to.  Can be the same as current.
    * @param genomesOnPath Intermediate genomes that must be visited
    * on the way down from coalescenceLimit to tgt.  If this is
    * specified as NULL, then the path will be computed automatically
    * (using hal::getGenomesInSpanningTree(coalescenceLimit, tgtGenome)).
    * Specifying this can avoid recomputing the path over and over again
    * when, say, calling getMappedSegments repeatedly for the same 
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
   virtual hal_size_t getMappedSegments(
     std::set<MappedSegmentConstPtr>& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath = NULL,
     bool doDupes = true,
     hal_size_t minLength = 0,
     const Genome *coalescenceLimit = NULL,
     const Genome *mrca = NULL) const = 0;

   /** Print contents of segment */
   virtual void print(std::ostream& os) const = 0;

protected:
   friend class counted_ptr<Segment>;
   friend class counted_ptr<const Segment>;
   virtual ~Segment() = 0;
};

inline Segment::~Segment() {}

inline std::ostream& operator<<(std::ostream& os, const Segment& s)
{
  s.print(os);
  return os;
}


}
#endif
