/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGENOME_H
#define _HALGENOME_H

#include <vector>
#include <string>
#include "halDefs.h"
#include "halSegmentedSequence.h"
#include "halSequence.h"

namespace hal {

/** 
 * Interface for a genome within a hal alignment.  The genome
 * is comprised of a dna sequence, and two segment arrays (top and bottom)
 * which are used to map between ancestral and descendant genomes.  This
 * data is all accessed through iterators. 
 */
class Genome : public SegmentedSequence
{
public:   
   /** Get the name of the genome */
   virtual const std::string& getName() const = 0;

   /** Reset (or initialize) the dimensions of the genome 
    * Note that there are no guarantees that any of the current
    * data gets preserved so this should only be used for creating
    * or completely rewriting a genome.  The phylogenetic information
    * (ex number of children) is read from the Aignment object
    * @param sequenceDimensions List of information for all sequences
    * that are to be contained in the genome. */
   virtual void setDimensions(
     const std::vector<hal::Sequence::Info>& sequenceDimensions) = 0;

   /** Reset (or initialize) the number of top segments in each
    * sequence of the genome, leaving the rest of the genome intact.
    * @param sequenceDimensions Number of top segments in each sequence
    * if there is a sequence presently in the genome but not in this 
    * list it will not be updated.*/
   virtual void setTopDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions) = 0;

   /** Reset (or initialize) the number of bottom segments in each
    * sequence of the genome, leaving the rest of the genome intact.
    * @param sequenceDimensions Number of bottom segments in each sequence
    * if there is a sequence presently in the genome but not in this 
    * list it will not be updated.*/
   virtual void setBottomDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions) = 0;
   
   /** Get number of sequences in the genome */
   virtual hal_size_t getNumSequences() const = 0;
   
   /** Get a sequence by name */
   virtual Sequence* getSequence(const std::string& name) = 0;

   /** Get a (read-only) sequence by name */
   virtual const Sequence* getSequence(const std::string& name) const = 0;

   /** Get a sequence by base's position (in genome coordinates) */
   virtual Sequence* getSequenceBySite(hal_size_t position) = 0;

   /** Get a sequence by base's position (in genome coordinates) */
   virtual const Sequence* getSequenceBySite(hal_size_t position) const = 0;
   
   /** Get a sequence iterator 
    * @param position Number of the sequence to start iterator at */
   virtual SequenceIteratorPtr getSequenceIterator(
     hal_index_t position = 0) = 0;

   /** Get a const sequence iterator 
    * @param position Number of the sequence to start iterator at */
   virtual SequenceIteratorConstPtr getSequenceIterator(
     hal_index_t position = 0) const = 0;

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
   virtual Genome* getChild(hal_size_t childIdx) = 0;

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
