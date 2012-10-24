/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLIFTOVER_H
#define _HALLIFTOVER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "hal.h"

namespace hal {

class Liftover
{
public:
   
   Liftover();
   virtual ~Liftover();

   void convert(AlignmentConstPtr alignment,
                const Genome* srcGenome,
                std::ifstream* inputFile,
                const Genome* tgtGenome,
                std::ofstream* outputFile);
                   
protected:

   bool readBedLine();
   void writeBedLine();
   void liftInterval();

   typedef ColumnIterator::DNASet DNASet;
   typedef ColumnIterator::ColumnMap ColumnMap;
   typedef PositionCache::IntervalSet IntervalSet;

protected: 

   AlignmentConstPtr _alignment;
   std::ifstream* _inputFile;
   std::ofstream* _outputFile;
   
   // current bed line
   std::string _inName;
   hal_size_t _inStart;
   hal_size_t _inEnd;
   std::vector<std::string> _data;

   std::string _buffer;
   std::string _outName;
   hal_size_t _outStart;
   hal_size_t _outEnd;

   const Genome* _srcGenome;
   const Genome* _tgtGenome;
   const Sequence* _srcSequence;
   std::set<const Genome*> _tgtSet;

   ColumnIteratorConstPtr _colIt;

   std::set<std::string> _missedSet;
};

}
#endif
