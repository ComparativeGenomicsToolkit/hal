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
   void convertBlock();
   void initArrayIndexes(size_t col);

   void aLine();
   void sLine();
   void end();

private:
   
   struct RowInfo 
   {
      hal_index_t _arrayIndex;
      size_t _gaps;
      const Record* _record;
   };

   std::string _refName;
   const DimMap* _dimMap;
   AlignmentPtr _alignment;
   std::vector<RowInfo> _blockInfo;
};

}

#endif
