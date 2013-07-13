#include "halPhyloP.h"

using namespace std;
using namespace hal;

halPhyloP::halPhyloP(hal::AlignmentConstPtr alignment, std::string modFile) {
  _softMaskDups = 1;  //use soft mask by default
  _maskAllDups = 0;  // only mask ambiguous bases in dups by default
  _alignment = alignment;
  _mode = CONACC;  // this tests for conservation and acceleration, I expect this is the only mode UCSC will want to use, but is easy to switch to CON or ACC or NNEUT

  //read in neutral model
  FILE *infile = phast_fopen(modFile.c_str(), "r");
  _mod = tm_new_from_file(infile, TRUE);
  phast_fclose(infile);

  //make sure all species in the tree are in the alignment, otherwise print warning and prune tree; create targetSet from these species. Make a hash of species names to species index (using phast's hash structure)
  int numspec=0;
  List *leafNames = tr_leaf_names(_mod->tree);
  int numleaf = lst_size(leafNames);
  List *pruneNames = lst_new_ptr(numleaf);
  _seqnameHash = hsh_new(numleaf*10);
  char **names = (char**)smalloc(numleaf * sizeof(char*));
  for (int i=0; i < lst_size(leafNames); i++) {
    string targetName = std::string(((String*)lst_get_ptr(leafNames, i))->chars);
    const Genome *tgtGenome = _alignment->openGenome(targetName);
    if (tgtGenome == NULL) {
      cerr << "Genome" << targetName << " not found in alignment; pruning from tree" << endl;
      lst_push(pruneNames, lst_get_ptr(leafNames, i));
    } else {
      String *leafName = (String*)lst_get_ptr(leafNames, i);
      hsh_put_int(_seqnameHash, leafName->chars, numspec);
      _targetSet.insert(tgtGenome);
      names[numspec] = (char*)smalloc((leafName->length+1)*sizeof(char));
      strcpy(names[numspec], leafName->chars);
      numspec++;
    }
  }
  lst_free_strings(leafNames);
  lst_free(leafNames);
  lst_free(pruneNames);

  //create a dummy alignment with a single column. As we iterate through the columns we will fill in the bases and compute the phyloP scores for individual columns.
  char **seqs = (char**)smalloc(numspec * sizeof(char*));
  for (int i=0; i < numspec; i++) {
    seqs[i] = (char*)smalloc(2*sizeof(char));
    seqs[i][0] = 'N';
    seqs[i][1] = '\0';
  }
  _msa = msa_new(seqs, names, numspec, 1, NULL);
  ss_from_msas(_msa, 1, 0, NULL, NULL, NULL, -1, FALSE);
  msa_free_seqs(_msa);
  _colfitdata = col_init_fit_data(_mod, _msa, ALL, _mode, FALSE);
  _colfitdata->tupleidx = 0;

}
				  
halPhyloP::~halPhyloP() {
  col_free_fit_data(_colfitdata);
  msa_free(_msa);
  hsh_free(_seqnameHash);
}

// set whether we want to use hard or soft masks for dups. Hard mask will turn the entire column to N's (and therefore p-value to 1) if any species is duplicated. Soft mask will turn just that species to N. I expect we only want soft mask, which is default.
void halPhyloP::setDupMask(string dupMask) {
  if (dupMask == "hard") {
    _softMaskDups=0;
  } else if (dupMask == "soft") {
    _softMaskDups=1;
  } else {
    throw hal_exception("unknown dupMask + " + dupMask + ", should be hard or soft");
  }
}

// set whether we want to consider anything a duplication if there are multiple alignments for one species, or if we only care about duplications which cause uncertainty in the bases (this is the default) 
void halPhyloP::setDupType(string dupType) {
  if (dupType == "ambiguous") {
    _maskAllDups = 0;
  } else if (dupType == "all") {
    _maskAllDups = 1;
  } else {
    throw hal_exception("unknown dupType + " + dupType + ", should be all or ambiguous");
  }
}

void halPhyloP::setPhyloPMode(string mode) {
  mode_type newmode;
  if (mode == "CONACC") {
    newmode = CONACC;
  } else if (mode == "CON") {
    newmode = CON;
  } else if (mode == "ACC") {
    newmode = ACC;
  } else if (mode == "NNEUT") {
    newmode = NNEUT;
  } else {
    throw hal_exception("unknown phyloP mode " + mode);
  }
  if (_mode != newmode) {
    col_free_fit_data(_colfitdata);
    _mode = newmode;
    _colfitdata = col_init_fit_data(_mod, _msa, ALL, _mode, FALSE);
  }
}


// compute phyloP score for a particular alignment column, return pval
double halPhyloP::pval(const ColumnIterator::ColumnMap *cmap) {
  for (int i=0; i < _msa->nseqs; i++) 
    _msa->ss->col_tuples[0][i]='*';
  
  for (ColumnIterator::ColumnMap::const_iterator it = cmap->begin(); it != cmap->end(); ++it) {
    const Sequence *sequence = it->first;
    const Genome *genome = sequence->getGenome();
    int spec = hsh_get_int(_seqnameHash, genome->getName().c_str());
    if (spec < 0)
      continue;
    ColumnIterator::DNASet* dnaSet = it->second;
    for (ColumnIterator::DNASet::const_iterator j = dnaSet->begin(); j != dnaSet->end(); ++j) {
      DNAIteratorConstPtr dna = *j;
      char base = toupper(dna->getChar());
      if (_msa->ss->col_tuples[0][spec]=='*')
	_msa->ss->col_tuples[0][spec] = base;
      else {
	if (_maskAllDups && _softMaskDups == 0) {  //hard mask, all dups
	  return 0.0;  // duplication; mask this base
	} else if (_maskAllDups) {  // soft mask, all dups
	  _msa->ss->col_tuples[0][spec] = 'N';
	} else if (_msa->ss->col_tuples[0][spec] != base) {
	  if (_softMaskDups == 0) 
	    return 0.0;
	  else _msa->ss->col_tuples[0][spec] ='N';
	} else {
	  _msa->ss->col_tuples[0][spec] = base;
	}
      }
    }
  }
  for (int i=0; i < _msa->nseqs; i++) {
    if (_msa->ss->col_tuples[0][i]=='*')
      _msa->ss->col_tuples[0][i]='N';
  }
  
  //finally, compute the score!
  double alt_lnl, null_lnl, this_scale, delta_lnl, pval;
  int sigfigs=4;  //same value used in phyloP code
  _mod->scale = 1;
  tm_set_subst_matrices(_mod);
  null_lnl = col_compute_log_likelihood(_mod, _msa, 0, _colfitdata->fels_scratch[0]);
  vec_set(_colfitdata->params, 0, _colfitdata->init_scale);
  opt_newton_1d(col_likelihood_wrapper_1d, &_colfitdata->params->data[0], _colfitdata,
		&alt_lnl, sigfigs, _colfitdata->lb->data[0], _colfitdata->ub->data[0],
		NULL, NULL, NULL);
  alt_lnl *= -1;
  this_scale = _colfitdata->params->data[0];
  delta_lnl = alt_lnl - null_lnl;
  if (delta_lnl <= -0.01) { //shouldn't happen
    cerr << "got delta_lnl = " << delta_lnl << endl;
    throw hal_exception(string("ERROR col_lrts: delta_lnl < 0 "));
  }
  if (_mode == NNEUT || _mode == CONACC) 
    pval = chisq_cdf(2*delta_lnl, 1, FALSE);
  else 
    pval = half_chisq_cdf(2*delta_lnl, 1, FALSE);
  /* assumes 50:50 mix of chisq and point mass at zero, due to
     bounding of param */
  
  if (pval < 1e-20)
    pval = 1e-20;
  /* approx limit of eval of tail prob; pvals of 0 cause problems */

  pval = -log10(pval);

  if (_mode == CONACC && this_scale > 1)
    pval *= -1; /* mark as acceleration */

  // print out parameters for debugging
//    cout << _msa->ss->col_tuples[0] << " " << _mod->scale << " " << alt_lnl << " " << null_lnl << " " << delta_lnl << " " << pval << " ";

  return pval;
}
    
