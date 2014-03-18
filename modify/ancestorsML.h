#ifndef __ANCESTORSML_H_
#define __ANCESTORSML_H_
#include "halDefs.h"
#include "halGenome.h"
#include "halAlignment.h"
#include "sonLibTree.h"

typedef struct {
  const hal::Genome *rootGenome;
  hal_index_t pos;
  // Reversed with respect to reference?
  bool reversed;
} rootInfo;

typedef struct {
  // Whether this node has already been calculated.
  bool done;
  // Position of this site in the genome
  hal_index_t pos;
  // Probability of leaves under this node given each nucleotide.
  double pLeaves[4];
  // points to maximum probability chars for the children of this node,
  // given each nucleotide for this node.
  char *MLChildrenChars[4];
  // array of probs for the chars of children of this node (given the assignment
  // of this node.)
  double *childrenCharProbs[4][4];
  // should only be set on the leaves at first.
  char dna;
  // Reversed with respect to reference?
  bool reversed;
  // phast ID from the model.
  int phastId;
} felsensteinData;

rootInfo * findRoot(const hal::Genome *genome,
                    hal_index_t pos, bool reversed=false);

void removeAndPrune(stTree *tree);

void buildTree(hal::AlignmentConstPtr alignment, const hal::Genome *genome,
               hal_index_t pos, stTree *tree, bool reversed, std::map<std::string, int> *nameToId = NULL);

void reEstimate(TreeModel *mod, hal::AlignmentPtr alignment, hal::Genome *genome, hal_index_t startPos, hal_index_t endPos, std::map<std::string, int> &nameToId, double threshold, bool writeHal, bool printWrites);

#endif
