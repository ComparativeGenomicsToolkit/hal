/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALWIGGLELIFTOVER_H
#define _HALWIGGLELIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "halWiggleScanner.h"
#include "halWiggleTiles.h"

namespace hal {

class WiggleLiftover : public WiggleScanner
{
public:
   
   WiggleLiftover();
   virtual ~WiggleLiftover();
   
   void preloadOutput(AlignmentConstPtr alignment,
                      const Genome* tgtGenome,
                      std::istream* inputFile);
   
   void convert(AlignmentConstPtr alignment,
                const Genome* srcGenome,
                std::istream* inputFile,
                const Genome* tgtGenome,
                std::ostream* outputFile,
                bool traverseDupes = true,
                bool unique = false);

   static const double DefaultValue;
   static const hal_size_t DefaultTileSize;

protected:

   virtual void visitLine();
   virtual void visitHeader();
   virtual void visitEOF();

   void mapSegment();
   void mapFragments(std::vector<MappedSegmentConstPtr>& fragments);
   void write();
                      
protected: 

   struct CoordVal { hal_index_t _first; hal_index_t _last; double _val; };
   typedef std::vector<CoordVal> ValVec;

   AlignmentConstPtr _alignment;
   std::istream* _inStream;
   std::ostream* _outStream;
   bool _traverseDupes;
   bool _unique;
      
   const Genome* _srcGenome;
   const Genome* _tgtGenome;
   const Sequence* _srcSequence;
   std::set<const Genome*> _tgtSet;
   std::set<MappedSegmentConstPtr> _mappedSegments;
   hal_index_t _lastIndex;

   SegmentIteratorConstPtr _segment;
   ValVec _cvals;
   WiggleTiles<double> _outVals;
   hal_index_t _cvIdx; 

};

}
#endif
