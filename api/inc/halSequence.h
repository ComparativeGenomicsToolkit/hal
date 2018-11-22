/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEQUENCE_H
#define _HALSEQUENCE_H

#include <string>
#include <set>
#include "halDefs.h"

namespace hal {

/** 
 * Interface for a sequence of DNA. Note that sequences should
 * not be written to outside of creating new genomes.
 */
class Sequence
{
public:
   /** Destructor */
    virtual ~Sequence() {
    }

   /** Basic information describing the dimenisons and name of a sequence
    * required for each sequence when creating a new genome */
   struct Info
   {
       Info():
        _length(0),
        _numTopSegments(0),
        _numBottomSegments(0) {}
      Info(const std::string& name,
           hal_size_t length,
           hal_size_t numTopSegments,
           hal_size_t numBottomSegments) :
        _name(name),
        _length(length),
        _numTopSegments(numTopSegments),
        _numBottomSegments(numBottomSegments) {}
      std::string _name;
      hal_size_t _length;
      hal_size_t _numTopSegments;
      hal_size_t _numBottomSegments;
   };

   /** Dimensional information for the number of (top or bottom) segments
    * in a sequence, required when partially updating a genome.  A partial
    * update in this context is a rewrite of either the bottom OR the top
    * segments and nothing else */
   struct UpdateInfo
   {
       UpdateInfo():
           _numSegments(0) {}
      UpdateInfo(const std::string& name,
                 hal_size_t numSegments) :
        _name(name),
        _numSegments(numSegments) {}
      std::string _name;
      hal_size_t _numSegments;
   };

   /** Return the sequence's name */
   virtual std::string getName() const = 0;

   /** Set the name of this sequence */
   virtual void setName(const std::string &newName) = 0;

   /** Return the sequence's name in genomeName.sequenceName format
    * (which is how it will appear in a MAF) */
   virtual std::string getFullName() const = 0;

    /** Get the total length of the DNA sequence in the sequence*/
   virtual hal_size_t getSequenceLength() const = 0;
   
   /** Get the containing (read-only) genome */
   virtual const Genome* getGenome() const = 0;

   /** Get the containing genome */
   virtual Genome* getGenome() = 0;

   /** Get the sequence's start position in the genome */
   virtual hal_index_t getStartPosition() const = 0;

   /** Get the sequence's end position (start + len - 1) in the genome */
   virtual hal_index_t getEndPosition() const = 0;

   /** Get the index of the sequence in the sequence array */
   virtual hal_index_t getArrayIndex() const = 0;

   /** Get the number of top segements 
    * (which form blocks with ancestor and siblings)
    * in the sequence */
   virtual hal_size_t getNumTopSegments() const = 0;

   /** Get the number of bottom segments
    * (which form blocks with the children)
    * in the sequence */
   virtual hal_size_t getNumBottomSegments() const = 0;

   /** Get the index of the sequence in the genome's top segment array */
   virtual hal_index_t getTopSegmentArrayIndex() const = 0;

   /** Get the index of the sequence in the genome's bottom segment array */
   virtual hal_index_t getBottomSegmentArrayIndex() const = 0;

   /** Get a top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position = 0) = 0;

   /** Get a const top segment iterator
    * @param position Index in segment array of returned iterator */
   virtual TopSegmentIteratorPtr getTopSegmentIterator(
     hal_index_t position = 0) const = 0;

   /** Get a bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position = 0) = 0;

   /** Get a const bottom segment iterator
    * @param position Index in segment array of returned iterator */
   virtual BottomSegmentIteratorPtr getBottomSegmentIterator(
     hal_index_t position = 0) const = 0;

   /** Create a DNA iterator
    * @param position Index in genome of returned iterator */
   virtual DnaIteratorPtr getDnaIterator(hal_index_t position = 0) = 0;

   /** Get a const DNA iterator
    * @param position Index in genome of returned iterator */
   virtual DnaIteratorPtr getDnaIterator(hal_index_t position = 0) const = 0;

   /** Get a column iterator that covers a specified range of this sequence.
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
   virtual ColumnIteratorPtr getColumnIterator(
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
     * sequence
     * @param outString String object into which we copy the result
     * @param start First position of substring 
     * @param length Length of substring */
    virtual void getSubString(std::string& outString, hal_size_t start,
                              hal_size_t length) const = 0;

  /** Set the character string underlying the sequence
    * @param inString input string to copy
    * @param start First position of substring 
    * @param length Length of substring */
   virtual void setSubString(const std::string& inString, 
                             hal_size_t start,
                             hal_size_t length) = 0;

   };

}
#endif
// Local Variables:
// mode: c++
// End:
