/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu) and 
 * Melissa Jane Hubisz (Cornell University)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halPhyloP.h"

using namespace std;
using namespace hal;

PhyloP::PhyloP() : _mod(NULL), _softMaskDups(false), _maskAllDups(false),
                   _seqnameHash(NULL), _colfitdata(NULL), _mode(CONACC),
                   _msa(NULL)
{
  
}

PhyloP::~PhyloP() 
{
  clear();
}

void PhyloP::clear()
{ 
  if (_colfitdata != NULL)
  {
    col_free_fit_data(_colfitdata);
  }
  if (_msa != NULL)
  {
    msa_free(_msa);
  }
  if (_seqnameHash != NULL)
  {
    hsh_free(_seqnameHash);
  }
  _targetSet.clear();

  // need to free _mod?

}

void PhyloP::init(AlignmentConstPtr alignment, const string& modFilePath,
                  ostream* outStream,
                  bool softMaskDups, 
                  const string& dupType,
                  const string& phyloPMode,
		  const string &subtree)
{
  clear();
  _alignment = alignment;
  _softMaskDups = (int)softMaskDups;
  _outStream = outStream;

  if (dupType == "ambiguous")
  {
    _maskAllDups = 0;
  }
  else if (dupType == "all")
  {
    _maskAllDups = 1;
  }
  else
  {
    throw hal_exception("unknown dupType + " + dupType + 
                        ", should be all or ambiguous");
  }

  if (phyloPMode == "CONACC")
  {
    _mode = CONACC;
  }
  else if (phyloPMode == "CON") 
  {
    _mode = CON;
  } 
  else if (phyloPMode == "ACC") 
  {
    _mode = ACC;
  }
  else if (phyloPMode == "NNEUT") 
  {
    _mode = NNEUT;
  }
  else 
  {
    throw hal_exception("unknown phyloP mode " + phyloPMode);
  }
  
  //read in neutral model
  FILE *infile = phast_fopen(modFilePath.c_str(), "r");
  _mod = tm_new_from_file(infile, TRUE);
  phast_fclose(infile);

 
  // make sure all species in the tree are in the alignment, otherwise print
  // warning and prune tree; create targetSet from these species.
  // Make a hash of species names to species index (using phast's hash
  // structure)
  int numspec=0;
  List *leafNames = tr_leaf_names(_mod->tree);
  int numleaf = lst_size(leafNames);
  List *pruneNames = lst_new_ptr(numleaf);
  _seqnameHash = hsh_new(numleaf*10);
  char **names = (char**)smalloc(numleaf * sizeof(char*));
  for (int i=0; i < lst_size(leafNames); i++) 
  {
    string targetName = string(((String*)lst_get_ptr(leafNames, i))->chars);
    const Genome *tgtGenome = _alignment->openGenome(targetName);
    if (tgtGenome == NULL) 
    {
      cerr << "Genome" << targetName 
           << " not found in alignment; pruning from tree" << endl;
      lst_push(pruneNames, lst_get_ptr(leafNames, i));
    } 
    else 
    {
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


  //create a dummy alignment with a single column. As we iterate through the 
  // columns we will fill in the bases and compute the phyloP scores for
  // individual columns.
  char **seqs = (char**)smalloc(numspec * sizeof(char*));
  for (int i = 0; i < numspec; i++) 
  {
    seqs[i] = (char*)smalloc(2*sizeof(char));
    seqs[i][0] = 'N';
    seqs[i][1] = '\0';
  }
  _msa = msa_new(seqs, names, numspec, 1, NULL);
  ss_from_msas(_msa, 1, 0, NULL, NULL, NULL, -1, FALSE);
  msa_free_seqs(_msa);


  if (subtree != "\"\"") {
    _mod->subtree_root = tr_get_node(_mod->tree, subtree.c_str());
    if (_mod->subtree_root == NULL) {
      tr_name_ancestors(_mod->tree);
      _mod->subtree_root = tr_get_node(_mod->tree, subtree.c_str());
      if (_mod->subtree_root == NULL)
	throw hal_exception("no node named " + subtree);
    }
    _modcpy = tm_create_copy(_mod);
    _modcpy->subtree_root = NULL;
    _colfitdata = col_init_fit_data(_modcpy, _msa, ALL, NNEUT, FALSE);
    _colfitdata2 = col_init_fit_data(_mod, _msa, SUBTREE, _mode, FALSE);
    _colfitdata2->tupleidx = 0;
    _insideNodes = lst_new_ptr(_mod->tree->nnodes);
    _outsideNodes = lst_new_ptr(_mod->tree->nnodes);
    tr_partition_leaves(_mod->tree, _mod->subtree_root, 
			_insideNodes, _outsideNodes);
  } else {
    _colfitdata = col_init_fit_data(_mod, _msa, ALL, _mode, FALSE);
  }
  _colfitdata->tupleidx = 0;
}

/** Given a Sequence (chromosome) and a (sequence-relative) coordinate
 * range, print the phyloP wiggle with respect to the genomes
 * in the target set 
 *
 * This is basically copied from halAlignability.cpp but returns
 * the phyloP score rather than the count of alignments.
 *
 * Coordinates are always genome-relative by default (as opposed to 
 * sequence-relative).  The one exception is all methods within the
 * Sequence interface. 
 *
 * By default, all bases in the reference genome are scanned.  And all
 * other genomes are considered.  The --refSequence, --start, and 
 * --length options can limit the query to a subrange.  Note that unless
 * --refSequence is specified, --start is genome-relative (based on 
 * all sequences being concatenated together).
 *
 * Genomes to include are determined by genomes in the neutral model 
 * (.mod) file- this should be in the format outputted by phyloFit.
 *
 * There are a few options for dealing with duplications. The default
 * is dupMask=soft, dupType=ambiguous. This replaces any species base 
 * with an N if duplications cause uncertainty in the base for a 
 * particular alignment column. If dupType=all, then all duplications
 * are masked regardless of whether they cause uncertainty. If
 * dupMask=hard, then the entire column is masked (p-value returned is 1.0)
 *
 * Print the phyloP wiggle for a subrange of a given sequence to
 * the output stream. */
void PhyloP::processSequence(const Sequence* sequence,
                             hal_index_t start,
                             hal_size_t length,
                             hal_size_t step)
{
  hal_size_t seqLen = sequence->getSequenceLength();
  if (seqLen == 0)
  {
    return;
  }
  /** If the length is 0, we do from the start position until the end
   * of the sequence */
  if (length == 0)
  {
    length = seqLen - start;
  }
  hal_size_t last = start + length;
  if (last > seqLen)
  {
    stringstream ss;
    ss << "Specified range [" << start << "," << length << "] is"
       << "out of range for sequence " << sequence->getName() 
       << ", which has length " << seqLen;
    throw (hal_exception(ss.str()));
  }

  const Genome* genome = sequence->getGenome();
  string sequenceName = sequence->getName();
  string genomeName = genome->getName();

  /** The ColumnIterator is fundamental structure used in this example to
   * traverse the alignment.  It essientially generates the multiple alignment
   * on the fly according to the given reference (in this case the target
   * sequence).  Since this is the sequence interface, the positions
   * are sequence relative.  Note that we must specify the last position
   * in advance when we get the iterator.  This will limit it following
   * duplications out of the desired range while we are iterating. */
  hal_size_t pos = start;
  ColumnIteratorConstPtr colIt = 
     sequence->getColumnIterator(&_targetSet,
                                 0, pos,
                                 last - 1);

  // note wig coordinates are 1-based for some reason so we shift to right
  *_outStream << "fixedStep chrom=" << sequenceName << " start=" << start + 1
              << " step=" << step << "\n";
  
  /** Since the column iterator stores coordinates in Genome coordinates
   * internally, we have to switch back to genome coordinates.  */
  // convert to genome coordinates
  pos += sequence->getStartPosition();
  last += sequence->getStartPosition();
  while (pos <= last)
  {
    /** ColumnIterator::ColumnMap maps a Sequence to a list of bases
     * the bases in the map form the alignment column.  Some sequences
     * in the map can have no bases (for efficiency reasons) */ 
    const ColumnIterator::ColumnMap* cmap = colIt->getColumnMap();
    double pval = this->pval(cmap);

    *_outStream << pval << '\n';
    
    /** lastColumn checks if we are at the last column (inclusive)
     * in range.  So we need to check at end of iteration instead
     * of beginning (which would be more convenient).  Need to 
     * merge global fix from other branch */
    if (colIt->lastColumn() == true)
    {
      break;
    }

    pos += step;    
    if (step == 1)
    {
      /** Move the iterator one position to the right */
      colIt->toRight();
      
      /** This is some tuning code that will probably be hidden from 
       * the interface at some point.  It is a good idea to use for now
       * though */
      // erase empty entries from the column.  helps when there are 
      // millions of sequences (ie from fastas with lots of scaffolds)
      if (pos % 1000 == 0)
      {
        colIt->defragment();
      }
    }
    else
    {
      /** Reset the iterator to a non-contiguous position */
      colIt->toSite(pos, last - 1);
    }
  }
}

// compute phyloP score for a particular alignment column, return pval
double PhyloP::pval(const ColumnIterator::ColumnMap *cmap) 
{
  for (int i=0; i < _msa->nseqs; i++) 
  {
    _msa->ss->col_tuples[0][i] = '*';
  }
  
  for (ColumnIterator::ColumnMap::const_iterator it = cmap->begin(); 
       it != cmap->end(); ++it) 
  {
    const Sequence *sequence = it->first;
    const Genome *genome = sequence->getGenome();
    int spec = hsh_get_int(_seqnameHash, genome->getName().c_str());
    if (spec < 0)
    {
      continue;
    }
    ColumnIterator::DNASet* dnaSet = it->second;
    for (ColumnIterator::DNASet::const_iterator j = dnaSet->begin(); 
         j != dnaSet->end(); ++j) 
    {
      DNAIteratorConstPtr dna = *j;
      char base = toupper(dna->getChar());
      if (_msa->ss->col_tuples[0][spec] == '*')
      {
	_msa->ss->col_tuples[0][spec] = base;
      }
      else 
      {
	if (_maskAllDups && _softMaskDups == 0) 
        {  //hard mask, all dups
	  return 0.0;  // duplication; mask this base
	} 
        else if (_maskAllDups) 
        {  // soft mask, all dups
	  _msa->ss->col_tuples[0][spec] = 'N';
	} 
        else if (_msa->ss->col_tuples[0][spec] != base) 
        {
	  if (_softMaskDups == 0) 
          {
	    return 0.0;
          }
	  else 
          {
            _msa->ss->col_tuples[0][spec] ='N';
          }
	} 
        else 
        {
	  _msa->ss->col_tuples[0][spec] = base;
	}
      }
    }
  }
  
  for (int i = 0; i < _msa->nseqs; i++) 
  {
    if (_msa->ss->col_tuples[0][i] == '*')
    {
      _msa->ss->col_tuples[0][i] = 'N';
    }
  }
  
  //finally, compute the score!
  double alt_lnl, null_lnl, this_scale, delta_lnl, pval;
  int sigfigs=4;  //same value used in phyloP code

  if (_mod->subtree_root == NULL) {
    _mod->scale = 1;
    tm_set_subst_matrices(_mod);
    null_lnl = col_compute_log_likelihood(_mod, _msa, 0,
					  _colfitdata->fels_scratch[0]);
    vec_set(_colfitdata->params, 0, _colfitdata->init_scale);
    opt_newton_1d(col_likelihood_wrapper_1d, &_colfitdata->params->data[0],
		  _colfitdata,
		  &alt_lnl, sigfigs, _colfitdata->lb->data[0], 
		  _colfitdata->ub->data[0],
		  NULL, NULL, NULL);
    alt_lnl *= -1;
    this_scale = _colfitdata->params->data[0];
  } else {  //subtree case
    if (!col_has_data_sub(_mod, _msa, 0, _insideNodes, _outsideNodes)) {
      alt_lnl = 0;
      null_lnl = 0;
      this_scale = 1;
    } else {
      vec_set(_colfitdata->params, 0, _colfitdata->init_scale);
      opt_newton_1d(col_likelihood_wrapper_1d, &_colfitdata->params->data[0],
		    _colfitdata,
		    &null_lnl, sigfigs, _colfitdata->lb->data[0], 
		    _colfitdata->ub->data[0], NULL, NULL, NULL);
      null_lnl *= -1;
      
      vec_set(_colfitdata2->params, 0, 
	      max(0.05, _colfitdata->params->data[0]));
      vec_set(_colfitdata2->params, 1, _colfitdata2->init_scale_sub);
      opt_bfgs(col_likelihood_wrapper, _colfitdata2->params, 
	       _colfitdata2, &alt_lnl, _colfitdata2->lb,
	       _colfitdata2->ub, NULL, NULL, OPT_HIGH_PREC, NULL, NULL);
      alt_lnl *= -1;
      this_scale = _colfitdata2->params->data[1];
    }
  }
  delta_lnl = alt_lnl - null_lnl;
  
  if (delta_lnl <= -0.01) 
  { 
    //shouldn't happen
    cerr << "got delta_lnl = " << delta_lnl << endl;
    throw hal_exception(string("ERROR col_lrts: delta_lnl < 0 "));
  }
  if (delta_lnl < 0) delta_lnl=0;
  if (_mode == NNEUT || _mode == CONACC) 
  {
    pval = chisq_cdf(2*delta_lnl, 1, FALSE);
  }
  else 
  {
    pval = half_chisq_cdf(2*delta_lnl, 1, FALSE);
  }
  /* assumes 50:50 mix of chisq and point mass at zero, due to
     bounding of param */
  
  if (pval < 1e-20)
  {
    pval = 1e-20;
  }
  /* approx limit of eval of tail prob; pvals of 0 cause problems */

  pval = -log10(pval);

  if (_mode == CONACC && this_scale > 1)
  {
    pval *= -1; /* mark as acceleration */
  }

  // print out parameters for debugging
//    cout << _msa->ss->col_tuples[0] << " " << _mod->scale << " " << alt_lnl << " " << null_lnl << " " << delta_lnl << " " << pval << " ";

  return pval;
}
    
