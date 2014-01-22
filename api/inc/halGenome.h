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
    * that are to be contained in the genome. 
    * @param storeDNAArrays Allocate arrays for DNA sequences.  This is
    * and should be true by default, but we introduce the option to allow
    * creating of HAL alignments without DNA information (ie just blocks). 
    * This functionality is used for halLodExtract, for example. */
   virtual void setDimensions(
     const std::vector<hal::Sequence::Info>& sequenceDimensions,
     bool storeDNAArrays = true) = 0;

   /** Update the number of top segments in *existing*
    * sequences of the genome, leaving the rest of the genome intact.
    * @param sequenceDimensions List of sequences names and their associated
    * number of top segments.  These sequences must have already been 
    * added using setDimensions (or error is thrown).*/
   virtual void updateTopDimensions(
     const std::vector<hal::Sequence::UpdateInfo>& sequenceDimensions) = 0;

   /** Update the number of bottom segments in *existing*
    * sequences of the genome, leaving the rest of the genome intact.
    * @param sequenceDimensions List of sequences names and their associated
    * number of bottom segments.  These sequences must have already been 
    * added using setDimensions (or error is thrown).*/
   virtual void updateBottomDimensions(
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

   /** Get a sequence end iterator (one beyond last element in list) 
    */
   virtual SequenceIteratorConstPtr getSequenceEndIterator() const = 0;

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

   /** Get the numeric index of a given child Genome 
    * @child child genome */
   virtual hal_index_t getChildIndex(const Genome* child) const = 0;

   /** Test if the genome stores DNA sequence.  Will be true unless
    * storeDNAArrays was set to false in setDimensions */
   virtual bool containsDNAArray() const = 0;

   /** Get a pointer to the alignment object that contains the genome.
    * Be careful not to free this pointer or put it inside an 
    * AlignmentConstPtr object since its memory is already spoken for */
   virtual const Alignment* getAlignment() const = 0;

   /** Copy all information from this genome to another. The genomes
    * must be in different alignments. The genome must not have
    * uninitialized data.
    * @param dest Genome to be copied to */
   void copy(Genome *dest) const;

   /** Copy all dimensions (sequence, top, and bottom) from this
    * genome to another.
    * @param dest Genome to be copied to */
   void copyDimensions(Genome *dest) const;

   /** Copy top dimensions from this genome to another (the genomes can be in
    * different alignment)
    * @param dest Genome to be copied to */
   void copyTopDimensions(Genome *dest) const;

   /** Copy bottom dimensions from this genome to another (the genomes
    * can be in different alignments). Only the segments corresponding
    * to child genomes that share the same name in both alignments are
    * copied. The segments for any other children of the destination
    * genome will need to be filled in later to avoid an invalid hal
    * file.
    * @param dest Genome to be copied to */
   void copyBottomDimensions(Genome *dest) const;

   /** Copy top segments from this genome to another (the genomes can be in
    * different alignments)
    * @param dest Genome to be copied to */
   void copyTopSegments(Genome *dest) const;

   /** Copy bottom segments from this genome to another. The genomes
    * must be in different alignments. The child indices do not have
    * to be consistent between the genomes.
    * @param dest Genome to be copied to */
   void copyBottomSegments(Genome *dest) const;

   /** Copy sequence from this genome to another (the genomes can be
    * in different alignments).
    * @param dest Genome to be copied to */
   void copySequence(Genome *dest) const;

   /** Copy metadata from this genome to another (the genomes can be
    * in different alignments).
    * @param dest Genome to be copied to */
   void copyMetadata(Genome *dest) const;

   /** Recompute parse info for this genome. */
   void fixParseInfo();

protected:

   /** Destructor */
   virtual ~Genome() = 0;
};

inline Genome::~Genome() {}

}
#endif
