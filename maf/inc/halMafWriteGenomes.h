/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFWRITEGENOMES_H
#define _HALMAFWRITEGENOMES_H

#include <fstream>
#include <string>
#include <deque>
#include <cstdlib>
#include <vector>
#include <string>
#include "hal.h"
#include "halMafScanner.h"
#include "halMafScanDimensions.h"

namespace hal {

/** update the HAL graph from the dimension information scanned
 * from the input MAF file, thumbing our nose at the haters as we do it!
 * remember that we convert lines to forward coordates as we read them
 * (but leave the strand character unchanged as a reminder */
class MafWriteGenomes : private MafScanner
{
public:
   MafWriteGenomes();
   ~MafWriteGenomes();
   
   typedef MafScanDimensions::DimMap DimMap;
   typedef MafScanDimensions::Record Record;
   typedef MafScanDimensions::StartMap StartMap;
   typedef MafScanDimensions::FilePosition FilePosition;
   typedef MafScanDimensions::PosSet PosSet;
   typedef std::pair<DimMap::const_iterator, DimMap::const_iterator> MapRange;

   void convert(const std::string& mafPath,
                const std::string& refGenomeName,
                const std::set<std::string>& targets,
                const DimMap& dimMap,
                AlignmentPtr alignment);
                         
private:
   
   MapRange getRefSequences() const;
   MapRange getNextSequences(DimMap::const_iterator jprev) const;

   void createGenomes();
   void initGenomes();
   void convertBlock();
   void initBlockInfo(size_t col);
   void initParaMap();
   void convertSegments(size_t col);
   void updateParalogy(size_t i);
   void initEmptySegments();
   void updateRefParseInfo();

   void aLine();
   void sLine();
   void end();

private:
   
   struct RowInfo 
   {
      hal_index_t _arrayIndex;
      size_t _gaps;
      hal_size_t _start;
      hal_size_t _length;
      const Record* _record;
      Genome* _genome;
      bool _skip;
   };

   struct Paralogy
   {
      hal_index_t _start;
      size_t _row;
      bool operator<(const Paralogy& p) const { return _start < p._start; }
   };
   typedef std::set<Paralogy> ParaSet;
   typedef std::map<Genome*, ParaSet> ParaMap;

private:

   ParaSet::iterator circularNext(size_t row, ParaSet& paraSet, 
                                  ParaSet::iterator i);
   ParaSet::iterator circularPrev(size_t row, ParaSet& paraSet,
                                  ParaSet::iterator i);

private:

   std::string _refName;
   Genome* _refGenome;
   hal_index_t _refRow;
   const DimMap* _dimMap;
   AlignmentPtr _alignment;
   std::vector<RowInfo> _blockInfo;
   TopSegmentIteratorPtr _topSegment, _paraTop;
   BottomSegmentIteratorPtr _bottomSegment, _refBottom;
   std::map<Genome*, hal_size_t> _childIdxMap;
   ParaMap _paraMap;
};

}

#endif
