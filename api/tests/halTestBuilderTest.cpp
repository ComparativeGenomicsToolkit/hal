/*
 * Copyright (C) 2012-2020 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halTestBuilder.h"
#include "halApiTestSupport.h"
#include <iostream>

struct TestBuilderTestBasic : public AlignmentTest {
    void createCallBack(Alignment *alignment) {
        double branchLength = 1e-10;
        GenomeSpec grandpa("grandpa");
        SeqSpec gseq(grandpa, "gseq", 17);

        GenomeSpec dad("dad", "grandpa", branchLength);
        SeqSpec dseq(dad, "dseq", 18);

        GenomeSpec son1("son1", "dad", branchLength);
        SeqSpec s1seq(son1, "s1seq", 15);
        GenomeSpec son2("son2", "dad", branchLength);
        SeqSpec s2seq(son2, "s2seq", 15);

        GenomeSpecVec genomes = {grandpa, dad, son1, son2};
        SeqSpecVec seqs = {gseq, dseq, s1seq, s2seq};

        BlockSpec block1 = {
            RowSpec(gseq,  2, 12, '-', "GCTATCGGGG"),
            RowSpec(dseq,  0, 10, '+', "GCTATCGGGG"),
            RowSpec(s1seq, 0, 10, '+', "GCTATCGGGG"),
            RowSpec(s2seq, 0, 10, '+', "GCTATCGGGG")
        };
        BlockSpec block2 = {
            RowSpec(gseq,  12, 17, '-', "gcctc"),
            RowSpec(dseq,  11, 16, '+', "gcctc"),
            RowSpec(s1seq, 10, 15, '+', "gcctc"),
            RowSpec(s2seq, 10, 15, '+', "gcctc")
        };
        BlockSpecVec blocks = {block1, block2};
        
        TestBuilder builder;
        builder.add(genomes, seqs, blocks);
        builder.build(alignment);
        builder.print(cerr);
    }
    void checkCallBack(const Alignment *alignment) {
    }
};

static void halTestBuilderTestBasic(CuTest *testCase) {
    TestBuilderTestBasic tester;
    tester.check(testCase);
}

static CuSuite *halTestBuilderTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halTestBuilderTestBasic);
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halTestBuilderTestSuite());
}
