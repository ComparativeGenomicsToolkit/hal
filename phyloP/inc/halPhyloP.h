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

class halPhyloP
{
 public:
  halPhyloP(hal::AlignmentConstPtr alignment, std::string modFile);
  ~halPhyloP();
  void setDupMask(std::string dupType);  // "hard" or "soft" dup maks
  void setDupType(std::string dupThreshold); // "ambiguous" or "all" dups
  void setPhyloPMode(std::string phyloPMode);  // "CONACC", "CON", "ACC", "NNEUT" are choices, though I think we are mainly interested in CONACC (conservation/acceleration- negative p-values indicate acceleration)
  double pval(const ColumnIterator::ColumnMap *cmap);  // return phyloP score 
  std::set<const Genome*> _targetSet;
 private:
  hal::AlignmentConstPtr _alignment;
  TreeModel* _mod;
  int _softMaskDups;  // 1 default = soft mask, if 0 use hard mask (mask entire column)
  int _maskAllDups;  // 0 default = mask only ambiguous bases in dups; if 1 mask any duplication
  hash_table* _seqnameHash;
  ColFitData* _colfitdata;
  mode_type _mode;
  MSA* _msa;
};
}
#endif
