/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halApiTestSupport.h"
#include "halAlignment.h"
#include "halGenome.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <unistd.h>
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

class AlignmentTestTrees : public AlignmentTest {
  public:
    void createCallBack(AlignmentPtr alignment) {
        hal_size_t alignmentSize = alignment->getNumGenomes();
        CuAssertTrue(_testCase, alignmentSize == 0);

        alignment->addRootGenome("Root", 0);
        alignment->addLeafGenome("Leaf", "Root", 10);
        alignment->addRootGenome("NewRoot", 15);
        alignment->addLeafGenome("Leaf1", "Root", 4.1);
        alignment->addLeafGenome("Leaf2", "Root", 5.1);
        alignment->addLeafGenome("Leaf3", "Root", 6.1);
        alignment->addLeafGenome("Leaf4", "Root", 7.1);
        alignment->updateBranchLength("Root", "Leaf1", 3.0);
        alignment->updateBranchLength("Root", "Leaf2", 6.1);
        alignment->updateBranchLength("Root", "Leaf2", 5.1);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        CuAssertTrue(_testCase, alignment->getRootName() == "NewRoot");
        CuAssertTrue(_testCase,
                     alignment->getNewickTree() == "((Leaf:10,Leaf1:3,Leaf2:5.1,Leaf3:6.1,Leaf4:7.1)Root:15)NewRoot;");
        CuAssertTrue(_testCase, alignment->getBranchLength("Root", "Leaf") == 10.0);
        CuAssertTrue(_testCase, alignment->getBranchLength("Root", "Leaf1") == 3.0);
        vector<string> children = alignment->getChildNames("Root");
        CuAssertTrue(_testCase, children.size() == 5);

        vector<string> leaves = alignment->getLeafNamesBelow("Leaf");
        CuAssertTrue(_testCase, leaves.size() == 0);
        leaves = alignment->getLeafNamesBelow("NewRoot");
        CuAssertTrue(_testCase, leaves.size() == 5);
        for (size_t i = 0; i < leaves.size(); ++i) {
            CuAssertTrue(_testCase, leaves[i][0] == 'L');
        }
        leaves = alignment->getLeafNamesBelow("Root");
        CuAssertTrue(_testCase, leaves.size() == 5);
        for (size_t i = 0; i < leaves.size(); ++i) {
            CuAssertTrue(_testCase, leaves[i][0] == 'L');
        }
    }
};

static void halAlignmentTestTrees(CuTest *testCase) {
    AlignmentTestTrees tester;
    tester.check(testCase);
}

static CuSuite *halAlignmentTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halAlignmentTestTrees);
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halAlignmentTestSuite());
}
