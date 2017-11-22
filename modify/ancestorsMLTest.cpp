#include <cstdio>
#include "hal.h"
#include "ancestorsML.h"

extern "C" {
#include "CuTest.h"
}

using namespace std;
using namespace hal;

// set up fake, test felsensteinData for this node.
static void labelTestTreeNode(stTree *node, map<string, int> *nameToId, char dna)
{
  felsensteinData *data = (felsensteinData *) malloc(sizeof(felsensteinData));
  memset(data, 0, sizeof(felsensteinData));
  if (nameToId != NULL) {
    data->phastId = (*nameToId)[string(stTree_getLabel(node))];
  }
  data->dna = dna;
  stTree_setClientData(node, data);
}

static void doFelsensteinWorkedExampleTest(CuTest *testCase)
{
  // Set up the PHAST substitution model.
  FILE *modFile = fopen("testdata/mammals.mod", "r");
  if (modFile == NULL) {
    throw hal_exception("can't find testdata/mammals.mod");
  }
  TreeModel *mod = tm_new_from_file(modFile, true);
  // Map names to phast model IDs.
  map<string, int> nameToId;
  List *phastList = tr_postorder(mod->tree);
  for (int i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = (TreeNode*) lst_get_ptr(phastList, i);
    nameToId[n->name] = n->id;
  }
  lst_free(phastList);
  tm_set_subst_matrices(mod);

  // Generate and label test tree with data A,G,A,A,C
  stTree *tree = stTree_parseNewickString("(((rat:0.2,mouse:0.2)mr:0.1,human:0.1)e:0.1,(cow:0.1,pig:0.1)l:0.1)b;");
  labelTestTreeNode(stTree_findChild(tree, "rat"), &nameToId, 'A');
  labelTestTreeNode(stTree_findChild(tree, "mouse"), &nameToId, 'G');
  labelTestTreeNode(stTree_findChild(tree, "human"), &nameToId, 'A');
  labelTestTreeNode(stTree_findChild(tree, "cow"), &nameToId, 'A');
  labelTestTreeNode(stTree_findChild(tree, "pig"), &nameToId, 'C');
  // Label ancestors with "unlabeled" base (Z)
  labelTestTreeNode(stTree_findChild(tree, "mr"), &nameToId, 'Z');
  labelTestTreeNode(stTree_findChild(tree, "e"), &nameToId, 'Z');
  labelTestTreeNode(stTree_findChild(tree, "l"), &nameToId, 'Z');
  labelTestTreeNode(stTree_findChild(tree, "b"), &nameToId, 'Z');

  doFelsenstein(tree, mod);
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  double likelihood = vec_get(mod->backgd_freqs, 0) * exp(data->pLeaves[0])
    // swap of 2, 1 here is intentional because 1 = G for us and 1 = C for the PHAST model
    + vec_get(mod->backgd_freqs, 2) * exp(data->pLeaves[1])
    + vec_get(mod->backgd_freqs, 1) * exp(data->pLeaves[2])
    + vec_get(mod->backgd_freqs, 3) * exp(data->pLeaves[3]);
  // phast's answer is -12.253583 -- check that we are reasonably close
  CuAssertDblEquals(testCase, -12.253583, log2(likelihood), 0.001);
}

int main(int argc, char *argv[]) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, doFelsensteinWorkedExampleTest);
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
  return suite->failCount > 0;
}
