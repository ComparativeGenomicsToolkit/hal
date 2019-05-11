/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALREARRANGEMENT_H
#define _HALREARRANGEMENT_H

#include <vector>
#include <string>
#include "halDefs.h"
#include "halGappedTopSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"

namespace hal {

/** 
 * Interface for a rearrangement operation.  Rearrangements are specifed by
 * pairs of breakpoints in a top segment array.  For now they are two-break
 * (balance or not) rearrangements and translocations.  Anything else will
 * will be labled complex for now. 
 */
class Rearrangement 
{
public:

   enum ID { Insertion = 0, Deletion, Duplication, Translocation, 
             Transposition, Inversion, Complex, Invalid, Gap, 
             Missing, Nothing};

   static const hal_size_t DefaultGapThreshold;
   static const double DefaultNThreshold;

   /* constructor */
   Rearrangement(const Genome* childGenome, 
                 hal_size_t gapThreshold = DefaultGapThreshold,
                 double nThreshold = DefaultNThreshold,
                 bool atomic = false);

   /** Destructor */
    virtual ~Rearrangement() {
    }

    /** Get the ID of the rearrangement */
   virtual ID getID() const;
   
   /** Get the length of a rearrangement (ie number of bases that were
    * inserted, deleted, moved, etc */
   virtual hal_size_t getLength() const;
   
   /** Get the number of gaps within the arrangement.  Gaps are dependent
    * on the gap threshold that can be specified below... */
   virtual hal_size_t getNumContainedGaps() const;

   /** Get the total size (in bases) of gaps within the arrangement.  
    * Gaps are dependent
    * on the gap threshold that can be specified below...*/
   virtual hal_size_t getNumContainedGapBases() const;

   /** Left breakpoint is specified by the first base of the returned segment */
   virtual TopSegmentIteratorPtr getLeftBreakpoint() const;

   /** Right breakpoint is specified by the last base of the returned
    * segment. If the rearrangement is a deletion, only the left
    * breakpoint is set. */
   virtual TopSegmentIteratorPtr getRightBreakpoint() const;

   /** Identify the rearrangement by scanning along the genome, skipping 
    * gaps until the matching breeakpoint is found.  If no matching breakpint
    * is found, the rearrangement cannot be identified and false is returned */
   virtual bool identifyFromLeftBreakpoint(TopSegmentIteratorPtr
                                           topSegment);

   /** Test if segment at given index corresponds to a deletion 
    * (shortcut method used by column iterator) */
   virtual bool identifyDeletionFromLeftBreakpoint(TopSegmentIteratorPtr
                                                   topSegment);

   /** Get the range in the parent that was deleted.  Only valid to call 
       this immediately after identifying a deletion with either of the
       above two functions.  Otherwise, results are undefined */
   virtual std::pair<hal_index_t, hal_index_t> getDeletedRange() const;

   /** Test if segment at given index corresponds to a breakpoint right
    * before an insertion 
    * (shortcut method used by column iterator) */
   virtual bool identifyInsertionFromLeftBreakpoint(TopSegmentIteratorPtr
                                                    topSegment);

   /** Get the range in the parent that was inserted.  Only valid to call 
    * this immediately after identifying a insertion.
    * Otherwise, results are undefined */
   virtual std::pair<hal_index_t, hal_index_t> getInsertedRange() const;

   /** Get the range in the parent that was duplicated.  Only valid to call 
    * this immediately after identifying a duplication.
    * Otherwise, results are undefined */
   virtual std::pair<hal_index_t, hal_index_t> getDuplicatedRange() const;

   /** Start scanning from last identified breakpoint */
   virtual bool identifyNext();

   /** Return the maximum size of a simple indel such that it will be 
    * considered a gap (and factored out of breakpoint computations) */
   virtual hal_size_t getGapLengthThreshold() const;

   /** Specify the maximum size of a simple indel such that it will be 
    * considered a gap (and factored out of breakpoint computations)  */
   virtual void setGapLengthThreshold(hal_size_t threshold);

   /** Set atomic behaviour (so gaps and compatible segments are not
    * used to merge larger segments.  Will overcount rearrangements but
    * useful sometimes (ie column iteraor).  Will automatically
    * set the gap threshold to 0.*/
   virtual void setAtomic(bool atomic);
   
   /* Get the atomic behaviour */
   virtual bool getAtomic() const;

   /* Set the threshold for missing data (N's) in a rearranged segment
    * @param nThreshold Value between 0 and 1.  If the number of N's
    * in the rearranged segment / segment length > nThreshold, then 
    * the rearrangement is termed "missing data".  */
   virtual void setNThreshold(double nThreshold);

   /* Get the threshold for missing data (N's) in a rearranged segment
    * If the number of N's
    * in the rearranged segment / segment length > nThreshold, then 
    * the rearrangement is termed "missing data". */
   virtual double getNThreshold() const;
   
private:
   
   virtual void resetStatus(TopSegmentIteratorPtr topSegment);
   
   virtual bool scanNothingCycle(TopSegmentIteratorPtr topSegment);
   virtual bool scanInversionCycle(TopSegmentIteratorPtr topSegment);
   virtual bool scanInsertionCycle(TopSegmentIteratorPtr topSegment);
   virtual bool scanDeletionCycle(TopSegmentIteratorPtr topSegment);
   virtual bool scanTranslocationCycle(TopSegmentIteratorPtr topSegment);
   virtual bool scanDuplicationCycle(TopSegmentIteratorPtr topSegment);

private:
   
   hal_size_t _gapThreshold;
   bool _atomic;
   double _nThreshold;

   GappedTopSegmentIteratorPtr _left, _cur, _next, _right;
   GappedBottomSegmentIteratorPtr _leftParent, _curParent, _nextParent,
   _rightParent;
   TopSegmentIteratorPtr _top;
   
   const Genome* _genome;
   const Genome* _parent;

   ID _id;  
};
}

#endif
// Local Variables:
// mode: c++
// End:
