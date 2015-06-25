/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTEDSEQUENCE_H
#define _HALSEGMENTEDSEQUENCE_H

#include "halDefs.h"
#include <set>

namespace hal {

/** 
 * Interface for a sequence of DNA that is also broken up into 
 * top and bottom segments.  This interface is extended by both
 * the Genome and Sequence classes. 
 *
 * Todo: Implenent the "Last" or maybe (right / left) iterator
 * notion, which will make looping easier and more general. 
 */
class SegmentedSequence
{
public:   

   /** Get the total length of the DNA sequence in the sequence*/
   virtual hal_size_t getSequenceLength() const = 0;
   
   /** Get the number of top segements 
    * (which form blocks with ancestor and siblings)
    * in the sequence */
   virtual hal_size_t getNumTopSegments() const = 0;

   /** Get the number of bottom segments
    * (which form blocks with the children)
    * in the sequence */
   virtual hal_size_t getNumBottomSegments() const = 0;

   /** Get a top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position = 0) = 0;

   /** Get a const top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorConstPtr getTopSegmentIterator(
     hal_index_t position = 0) const = 0;

   /** Get a topSegment end iterator (one beyond last element in list) */
   virtual TopSegmentIteratorConstPtr getTopSegmentEndIterator() const = 0;

   /** Get a bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position = 0) = 0;

   /** Get a const bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorConstPtr getBottomSegmentIterator(
     hal_index_t position = 0) const = 0;

   /** Get a bottomSegment end iterator (one beyond last element in list) */
   virtual BottomSegmentIteratorConstPtr getBottomSegmentEndIterator() 
     const = 0;

   /** Get a DNA iterator
    * @param position Index in genome of returned iterator */
   virtual DNAIteratorPtr getDNAIterator(
     hal_index_t position = 0) = 0;

   /** Get a const DNA iterator
    * @param position Index in genome of returned iterator */
   virtual DNAIteratorConstPtr getDNAIterator(
     hal_index_t position = 0) const = 0;

   /** Get a DNA end iterator (one beyond last element in list) */
   virtual DNAIteratorConstPtr getDNAEndIterator() const = 0;

   /** Get a column iterator 
    * @param targets Only genomes in this set are visited
    * * (note that other genomes in their spanning tree will be
    * * traversed as necessary but not reported)
    * @param maxInsertLength maximum insertion to be traversed.
    * Note: if setting > 0, it is recommended to set unique to true.
    * @param position Index in sequence of returned iterator 
    * @param noDupes Don't follow paralogy edges
    * @param noAncestors Don't report any non-leaf nodes in output
    * @param reverseStrand Map from reverse strand of this sequence
    * (but still in a left-to-right direction of forward strand) 
    * @param unique calls to toRight() will automatically skip 
    * over bases in the reference sequence that have already been
    * seen in an alignment column (ie homologous to already visited
    * sites) using the visitCache.
    * @param onlyOrthologs Include only the orthologs for each
    * reference base. In practice, this means paralogy edges are only
    * followed when moving down the tree for a particular column,
    * never up. */
   virtual ColumnIteratorConstPtr getColumnIterator(
     const std::set<const Genome*>* targets = NULL,
     hal_size_t maxInsertLength = 0,
     hal_index_t position = 0,
     hal_index_t lastPosition = NULL_INDEX,
     bool noDupes = false,
     bool noAncestors = false,
     bool reverseStrand = false,
     bool unique = false,
     bool onlyOrthologs = false) const = 0;

   /** Get the character string underlying the segmented sequence
    * @param outString String object into which we copy the result */
   virtual void getString(std::string& outString) const = 0;

  /** Set the character string underlying the segmented sequence
    * @param inString input string to copy */
   virtual void setString(const std::string& inString) = 0;

   /** Get the substring of character string underlying the 
    * segmented sequence
    * @param outString String object into which we copy the result
    * @param start First position of substring 
    * @param length Length of substring */
   virtual void getSubString(std::string& outString, hal_size_t start,
                             hal_size_t length) const = 0;

  /** Set the character string underlying the segmented sequence
    * @param inString input string to copy
    * @param start First position of substring 
    * @param length Length of substring */
   virtual void setSubString(const std::string& inString, 
                             hal_size_t start,
                             hal_size_t length) = 0;

   /** Get a rearrangement object 
    * @param position Position of topsegment defining first breakpoint of
    * rearrangement 
    * @param gapLengthThreshold The maximum size of a simple indel such that
    * it will be considered a gap (and factored out of breakpoint 
    * computations)
    * @param nThreshold Maximum fraction (between 0 and 1) of N's in a 
    * rearrangement for it not to be flagged as missing data
    * @param atomic Never merge segments together into a rearrangement.  
    * Generally for internal use */
   virtual RearrangementPtr getRearrangement(hal_index_t position,
                                             hal_size_t gapLengthThreshold,
                                             double nThreshold,
                                             bool atomic = false) const = 0;

   /** Get a gapped iterator  (experimental) 
    * REMINDER: ther iterator is extended from the left-to-right along
    * the sequence from i.  A seperate function is needed to, say,
    * get the last iterator from a sequence  (not a big deal since
    * the functinality is already in the implementation (extendLeft)*/
   virtual GappedTopSegmentIteratorConstPtr getGappedTopSegmentIterator(
     hal_index_t i, hal_size_t gapThreshold, bool atomic = false) const = 0;

   /** Get a gapped iterator  (experimental) 
    * REMINDER: ther iterator is extended from the left-to-right along
    * the sequence from i.  A seperate function is needed to, say,
    * get the last iterator from a sequence  (not a big deal since
    * the functinality is already in the implementation (extendLeft)*/
   virtual GappedBottomSegmentIteratorConstPtr getGappedBottomSegmentIterator(
     hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
     bool atomic = false) const = 0;

protected:

   /** Destructor */
   virtual ~SegmentedSequence() = 0;
};

inline SegmentedSequence::~SegmentedSequence() {}

}
#endif
