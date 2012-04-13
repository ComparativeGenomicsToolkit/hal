/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGENOME_H
#define _HALGENOME_H

#include "halDefs.h"

namespace hal {

/** 
 * Interface for a genome within a hal alignment.  The genome
 * is comprised of a dna sequence, and two segment arrays (top and bottom)
 * which are used to map between ancestral and descendant genomes.  This
 * data is all accessed through iterators. 
 */
class Genome
{
public:   

   /** Reset (or initialize) the dimensions of the genome 
    * Note that there are no guarantees that any of the current
    * data gets preserved so this should only be used for creating
    * or completely rewriting a genome.  The phylogenetic information
    * (ex number of children) is read from the Aignment object
    * @param totalSequenceLength DNA sequence length 
    * @param numTopSegments Number of top segments
    * @param numBottomSegments Number of bottom segments */
   virtual void reset(hal_size_t totalSequenceLength,
                      hal_size_t numTopSegments,
                      hal_size_t numBottomSegments) = 0;

   /** Reset (or initialize) the top segments, leaving the 
    * rest of the genome intact.
    * @param numTopSegments Number of top segments */
   virtual void resetTopSegments(hal_size_t numTopSegments) = 0;

   /** Reset (or initialize) the bottom segments, leaving the 
    * rest of the genome intact.
    * @param numBottomSegments Number of bottom segments */
   virtual void resetBottomSegments(hal_size_t numBottomSegments) = 0;
     
   /** Get the name of the genome */
   virtual const std::string& getName() const = 0;
   
   /** Get the total length of the DNA sequence in the genome*/
   virtual hal_size_t getSequenceLength() const = 0;
   
   /** Get the number of top segements 
    * (which form blocks with ancestor and siblings)
    * in the genome */
   virtual hal_size_t getNumberTopSegments() const = 0;

   /** Get the number of bottom segments
    * (which form blocks with the children)
    * in the genome */
   virtual hal_size_t getNumberBottomSegments() const = 0;

   /** Get a top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position) = 0;

   /** Get a const top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorConstPtr getTopSegmentIterator(
     hal_index_t position) const = 0;

   /** Get a bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position) = 0;

   /** Get a const bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorConstPtr getBottomSegmentIterator(
     hal_index_t position) const = 0;

   /** Get genome-specific metadata for this genome */
   virtual MetaData* getMetaData() = 0;

   /** Get read-only instance of genome-specific metadata for this genome */
   virtual const MetaData* getMetaData() const = 0;

   /** Get the parent genome */
   virtual Genome* getParent() = 0;

   /** Get const parent genome */
   virtual const Genome* getParent() const = 0;

   /** Get child genome 
    * @param childIdx index of child genome */
   virtual const Genome* getChild(hal_size_t childIdx) = 0;

   /** Get const child genome 
    * @param childIdx index of child genome */
   virtual const Genome* getChild(hal_size_t childIdx) const = 0;

   /** Get number of child genomes */
   virtual hal_size_t getNumChildren() const = 0;

protected:

   /** Destructor */
   virtual ~Genome() = 0;
};

inline Genome::~Genome() {}

}
#endif
