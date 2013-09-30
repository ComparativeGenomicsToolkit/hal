/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEQUENCE_H
#define _HALSEQUENCE_H

#include <string>
#include "halDefs.h"
#include "halSegmentedSequence.h"

namespace hal {

/** 
 * Interface for a sequence of DNA. Note that sequences should
 * not be written to outside of creating new genomes.
 */
class Sequence : public SegmentedSequence
{
public:

   /** Basic information describing the dimenisons and name of a sequence
    * required for each sequence when creating a new genome */
   struct Info
   {
      Info(){}
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
      UpdateInfo(){}
      UpdateInfo(const std::string& name,
                 hal_size_t numSegments) :
        _name(name),
        _numSegments(numSegments) {}
      std::string _name;
      hal_size_t _numSegments;
   };

   /** Return the sequence's name */
   virtual std::string getName() const = 0;

   /** Return the sequence's name in genomeName.sequenceName format
    * (which is how it will appear in a MAF) */
   virtual std::string getFullName() const = 0;

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

   /** Get the index of the sequence in the genome's top segment array */
   virtual hal_index_t getTopSegmentArrayIndex() const = 0;

   /** Get the index of the sequence in the genome's bottom segment array */
   virtual hal_index_t getBottomSegmentArrayIndex() const = 0;

protected:
   
   /** Destructor */
   virtual ~Sequence() = 0;
};

inline Sequence::~Sequence() {}

}
#endif
