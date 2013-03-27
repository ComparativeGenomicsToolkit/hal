/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <deque>
#include "halLodExtract.h"
extern "C" {
#include "sonLibTree.h"
}
 

using namespace std;
using namespace hal;

LodExtract::LodExtract()
{

}

LodExtract::~LodExtract()
{
  
}

void LodExtract::createInterpolatedAlignment(AlignmentConstPtr inAlignment,
                                             AlignmentPtr outAlignment,
                                             hal_size_t step,
                                             const string& tree)
{
  _inAlignment = inAlignment;
  _outAlignment = outAlignment;
  
  string newTree = tree.empty() ? inAlignment->getNewickTree() : tree;
  createTree(newTree);
}

void LodExtract::createTree(const string& tree)
{
  if (_outAlignment->getNumGenomes() != 0)
  {
    throw hal_exception("Output alignment not empty");
  }
  stTree* root = stTree_parseNewickString(const_cast<char*>(tree.c_str()));
  assert(root != NULL);
  deque<stTree*> bfQueue;
  bfQueue.push_front(root);
  while (!bfQueue.empty())
  {
    stTree* node = bfQueue.back();
    const char* label = stTree_getLabel(node);
    if (label == NULL)
    {
      throw hal_exception("Error parsing tree: unlabeled node");
    }
    const Genome* test = _inAlignment->openGenome(label);
    if (test == NULL)
    {
      throw hal_exception(string("Genome in tree: ") + string(label) +
                          "doesn't exist in source alignment");
    }
    else
    {
      _inAlignment->closeGenome(test);
    }    
    if (node == root)
    {
      _outAlignment->addRootGenome(label);
    }
    else
    {
      stTree* parent = stTree_getParent(node);
      assert(parent != NULL);
      const char* pLabel = stTree_getLabel(parent);
      assert(pLabel != NULL);
      double branchLength = stTree_getBranchLength(node);
      // clamp undefined branch lengths to 1. for now
      if (branchLength > 1e10)
      {
        branchLength = 1.;
      }
      _outAlignment->addLeafGenome(label, pLabel, branchLength);
    }

    int32_t numChildren = stTree_getChildNumber(node);
    for (int32_t childIdx = 0; childIdx < numChildren; ++childIdx)
    {
      bfQueue.push_front(stTree_getChild(node, childIdx));
    }

    bfQueue.pop_back();
    stTree_destruct(node);
  }
}
