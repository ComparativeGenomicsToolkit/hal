/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu) and 
 * Melissa Jane Hubisz (Cornell University)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALPHYLOP_H
#define _HALPHYLOP_H

#include <cstdlib>
#include <string>
#include "hal.h"

extern "C"{
#include "tree_model.h"
#include "fit_column.h"
#include "msa.h"
#include "sufficient_stats.h"
#include "hashtable.h"
}

namespace hal {

/** Use the Phast library methods to compute a PhyloP score for a HAL
 * aligment, column by column.  Thanks to Melissa Jane Hubisz. */
class PhyloP
{
public:

   PhyloP();
   virtual ~PhyloP();
   
   /** @param dupHardMask true for hard duplication mask or false for 
    * soft duplication mask
    * @param dupType ambiguous or all
    * @param phyloPMode "CONACC", "CON", "ACC", "NNEUT" are choices, 
    * though I think we are mainly interested in CONACC 
    * (conservation/acceleration- negative p-values indicate acceleration)
    * @param subtree If equal to empty string, perform phyloP test on 
    * entire tree. Otherwise, subtree names a branch to perform test on
    * subtree relative to rest of tree. The subtree includes all children
    * of the named node as well as the branch leading to the node.
    */
   void init(AlignmentConstPtr alignment, const std::string& modFilePath,
             std::ostream* outStream,
             bool softMaskDups = true, 
             const std::string& dupType = "ambiguous",
             const std::string& phyloPMode = "CONACC",
             const std::string& subtree = "\"\"");

   void processSequence(const Sequence* sequence,
                        hal_index_t start,
                        hal_size_t length,
                        hal_size_t step);

protected:

   // return phyloP score 
   double pval(const ColumnIterator::ColumnMap *cmap);  

   void clear();

protected:

   hal::AlignmentConstPtr _alignment;
   TreeModel* _mod;
   TreeModel* _modcpy;
   std::set<const Genome*> _targetSet;
   std::ostream* _outStream;
  
   // 1 default = soft mask, if 0 use hard mask (mask entire column)
   int _softMaskDups;  
   int _maskAllDups;  

   // 0 default = mask only ambiguous bases in dups; if 1 mask any duplication
   hash_table* _seqnameHash;
   ColFitData* _colfitdata;
   ColFitData *_colfitdata2;
   List *_insideNodes;
   List *_outsideNodes;
   mode_type _mode;
   MSA* _msa;
};

}
#endif
