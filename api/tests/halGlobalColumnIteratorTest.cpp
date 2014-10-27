#include "halAlignmentTest.h"
#include "halAlignment.h"
#include "halGlobalColumnIterator.h"
#include "halGlobalColumnIteratorTest.h"
#include "halColumnIterator.h"
#include "halRandomData.h"

using namespace std;
using namespace hal;

// Tests properties that should hold for all alignments on a randomly
// generated alignment.

void GlobalColumnIteratorBasicTest::createCallBack(AlignmentPtr alignment)
{
    createRandomAlignment(alignment, 10, 1e-10, 3, 77, 77, 10, 10);
}

void GlobalColumnIteratorBasicTest::checkCallBack(AlignmentConstPtr alignment)
{
    GlobalColumnIteratorPtr colIt = alignment->getGlobalColumnIterator(
        NULL, // all genomes
        false, // show dups
        false, // show ancestors
        false); // don't reverse strand
    // For ensuring that no column is repeated.
    ColumnIterator::VisitCache visitCache;

    for (;;) {
        const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
        for (ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
             colMapIt != colMap->end(); colMapIt++) {
            const Genome *genome = colMapIt->first->getGenome();
            for (hal_size_t i = 0; i < colMapIt->second->size(); i++) {
                DNAIteratorConstPtr dnaIt = colMapIt->second->at(i);
                hal_index_t pos = dnaIt->getArrayIndex();
                // Make sure we haven't seen this position before
                CuAssertTrue(_testCase, !visitCache[genome]->find(pos));
                // Record this position
                visitCache[genome]->insert(pos);
            }
        }
        if (colIt->lastColumn()) {
            break;
        }
        colIt->toRight();
    }

    // Make sure we've visited every position in the leaves of the alignment.
    vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
    vector<const Genome *> leafGenomes;
    for (hal_size_t i = 0; i < leafNames.size(); i++) {
        leafGenomes.push_back(alignment->openGenome(leafNames[i]));
    }
    for (hal_size_t i = 0; i < leafGenomes.size(); i++) {
        CuAssertTrue(_testCase, visitCache[leafGenomes[i]]->size() == leafGenomes[i]->getSequenceLength());
    }
}

void halGlobalColumnIteratorBasicTest(CuTest *testCase)
{
    try {
        GlobalColumnIteratorBasicTest test;
        test.check(testCase);
    } catch (...) {
        CuAssertTrue(testCase, false);
    }
}

CuSuite *halGlobalColumnIteratorTestSuite(void)
{
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halGlobalColumnIteratorBasicTest);
    return suite;
}
