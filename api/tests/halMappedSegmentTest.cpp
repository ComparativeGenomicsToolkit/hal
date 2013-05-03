/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "halMappedSegmentTest.h"
#include "halRandomData.h"
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"
#include "hal.h"


using namespace std;
using namespace hal;

void MappedSegmentMapUpTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // setup simple case were there is an edge from a parent to 
  // child and it is reversed
  Genome* parent = alignment->addRootGenome("parent");
  Genome* child1 = alignment->addLeafGenome("child1", "parent", 1);
  Genome* child2 = alignment->addLeafGenome("child2", "parent", 1);
  seqVec[0] = Sequence::Info("Sequence", 10, 0, 1);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 10, 1, 0);
  child1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 10, 1, 0);
  child2->setDimensions(seqVec);

  parent->setString("CCCTACGTGC");
  child1->setString("CCCTACGTGC");
  child2->setString("CCCTACGTGC");

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 10, 0);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 10, 0, true, 0);
  ts.applyTo(ti);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 10, 0, false, 0);
  ts.applyTo(ti);
}

void MappedSegmentMapUpTest::testTopSegment(AlignmentConstPtr alignment,
                                            TopSegmentIteratorConstPtr top)
{
  const Genome* parent = top->getGenome()->getParent();
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, parent, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  CuAssertTrue(_testCase, mseg->getGenome() == bottom->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == bottom->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == bottom->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == bottom->getReversed());
}

void MappedSegmentMapUpTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");
  TopSegmentIteratorConstPtr top = child2->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
  top = child1->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
}

void MappedSegmentMapDownTest::testBottomSegment(
  AlignmentConstPtr alignment,
  BottomSegmentIteratorConstPtr bottom,
  hal_size_t childIndex)
{
  const Genome* child = bottom->getGenome()->getChild(childIndex);
  vector<MappedSegmentConstPtr> results;
  bottom->getMappedSegments(results, child, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == 
               bottom->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               bottom->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == bottom->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == bottom->getReversed());
  TopSegmentIteratorConstPtr top = child->getTopSegmentIterator();
  top->toChild(bottom, childIndex);
  CuAssertTrue(_testCase, mseg->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == top->getReversed());
}

void MappedSegmentMapDownTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");

  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
  bottom->slice(1,2);
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
  bottom->toReverse();
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
}

void MappedSegmentMapAcrossTest::testTopSegment(AlignmentConstPtr alignment,
                                                TopSegmentIteratorConstPtr top)
{
  const Genome* parent = top->getGenome()->getParent();
  const Genome* other = top->getGenome()->getName() == "child1" ? 
     alignment->openGenome("child2") : alignment->openGenome("child1");
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, other, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  TopSegmentIteratorConstPtr sister = other->getTopSegmentIterator();
  sister->toChildG(bottom, other);
  CuAssertTrue(_testCase, mseg->getGenome() == sister->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == sister->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == sister->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == sister->getReversed());
}

void MappedSegmentMapAcrossTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");
  TopSegmentIteratorConstPtr top = child2->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
  top = child1->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
}

void MappedSegmentMapDupeTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // setup simple case were there is an edge from a parent to 
  // child and it is reversed
  Genome* parent = alignment->addRootGenome("parent");
  Genome* child1 = alignment->addLeafGenome("child1", "parent", 1);
  Genome* child2 = alignment->addLeafGenome("child2", "parent", 1);
  seqVec[0] = Sequence::Info("Sequence", 3, 0, 1);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  child1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  child2->setDimensions(seqVec);

  parent->setString("CCC");
  child1->setString("CCCTACGTG");
  child2->setString("CCCTACGTG");

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 3, 0);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 3, 0, true, 0, 1);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, true, 0, 2);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, 0, true, 0, 0);
  ts.applyTo(ti);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 3, 0, false, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, true, NULL_INDEX);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, 0, false, NULL_INDEX);
  ts.applyTo(ti);
}

void MappedSegmentMapDupeTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");

  TopSegmentIteratorConstPtr top = child1->getTopSegmentIterator();
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, child2, NULL, true);
  CuAssertTrue(_testCase, results.size() == 1);
  
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  TopSegmentIteratorConstPtr sister = child2->getTopSegmentIterator();
  sister->toChildG(bottom, child2);
  CuAssertTrue(_testCase, mseg->getGenome() == sister->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == sister->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == sister->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == sister->getReversed());

  top = child2->getTopSegmentIterator();
  results.clear();
  sister = child1->getTopSegmentIterator();
  top->getMappedSegments(results, child1, NULL, true);
  CuAssertTrue(_testCase, results.size() == 3);
  bool found[3] = {false};
  vector<MappedSegmentConstPtr>::iterator i = results.begin();
  for (; i != results.end(); ++i)
  {
    MappedSegmentConstPtr mseg = *i;
    CuAssertTrue(_testCase, mseg->getSource()->getGenome() == 
                 top->getGenome());
    CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
                 top->getStartPosition());
    CuAssertTrue(_testCase, 
                 mseg->getSource()->getLength() == top->getLength());
    CuAssertTrue(_testCase, 
                 mseg->getSource()->getReversed() == top->getReversed());
    BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
    bottom->toParent(top);
    TopSegmentIteratorConstPtr sister = child2->getTopSegmentIterator();
    sister->toChildG(bottom, child1);
    CuAssertTrue(_testCase, mseg->getGenome() == sister->getGenome());
    CuAssertTrue(_testCase, 
                 mseg->getLength() == sister->getLength());
    found[mseg->getArrayIndex()] = true;
  }
  CuAssertTrue(_testCase, found[0] == true);
  CuAssertTrue(_testCase, found[1] == true);
  CuAssertTrue(_testCase, found[2] == true);
}

void  MappedSegmentColCompareTest::checkCallBack(AlignmentConstPtr alignment)
{
  if (alignment->getNumGenomes() == 0)
  {
    return;
  }
  set<const Genome*> genomeSet;
  hal::getGenomesInSubTree(alignment->openGenome(alignment->getRootName()), 
                           genomeSet);
  for (set<const Genome*>::iterator i = genomeSet.begin(); i != genomeSet.end();
       ++i)
  {
    const Genome* srcGenome = *i;
    for (set<const Genome*>::iterator j = genomeSet.begin(); 
         j != genomeSet.end(); ++j)
    {
      const Genome* tgtGenome = *j;

      if (srcGenome->getSequenceLength() > 0 && 
          tgtGenome->getSequenceLength() > 0)
      {
        createColArray(srcGenome, tgtGenome);
        createBlockArray(srcGenome, tgtGenome);
        compareArrays();
      }
    }
  }
}

void  MappedSegmentColCompareTest::createColArray(const Genome* ref, 
                                                  const Genome* tgt)
{
  hal_size_t N = ref->getSequenceLength();
  _colArray.clear();
  _colArray.resize(N);
  set<const Genome*> tgtSet;
  tgtSet.insert(tgt);
  ColumnIteratorConstPtr colIt = ref->getColumnIterator(&tgtSet);
  for (; colIt->lastColumn() == false; colIt->toRight())
  {
    const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
    ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
    for (; colMapIt != colMap->end(); colMapIt++)
    {
      if (colMapIt->first->getGenome() == tgt)
      {
        ColumnIterator::DNASet* dnaSet = colMapIt->second;
        ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
        for (; dnaIt != dnaSet->end(); ++dnaIt)
        {
          DNAIteratorConstPtr dna = *dnaIt;
          _colArray[colIt->getReferenceSequencePosition()].insert(
            pair<hal_index_t, bool>(dna->getArrayIndex(), dna->getReversed()));
        }
      }
      else
      {
        CuAssertTrue(_testCase, colMapIt->first->getGenome() == ref);
      }
    }
  }
}

void  MappedSegmentColCompareTest::createBlockArray(const Genome* ref, 
                                                    const Genome* tgt)
{
  hal_size_t N = ref->getSequenceLength();
  _blockArray.clear();
  _blockArray.resize(N);
  SegmentIteratorConstPtr refSeg;
  hal_index_t numSegs;
  if (ref->getNumTopSegments() > 0)
  { 
    refSeg = ref->getTopSegmentIterator(0);
    numSegs = ref->getNumTopSegments();
  }
  else
  {
    refSeg = ref->getBottomSegmentIterator(0);
    numSegs = ref->getNumBottomSegments();
  }

  vector<MappedSegmentConstPtr> results;  
  for (; refSeg->getArrayIndex() < numSegs; refSeg->toRight())
  {
    refSeg->getMappedSegments(results, tgt);
  }
  
  for (vector<MappedSegmentConstPtr>::iterator i = results.begin();
       i != results.end(); ++i)
  {
    MappedSegmentConstPtr mseg = *i;
    SlicedSegmentConstPtr refSeg = mseg->getSource();
    hal_index_t refDelta = refSeg->getReversed() ? -1 : 1;
    hal_index_t mDelta = mseg->getReversed() ? -1 : 1;
    CuAssertTrue(_testCase, refSeg->getGenome() == ref);

    for (hal_index_t offset = 0; offset < mseg->getLength(); ++offset)
    {
      hal_index_t refPos = refSeg->getStartPosition() + offset * refDelta;
      hal_index_t mPos = mseg->getStartPosition() + offset * mDelta;
      _blockArray[refPos].insert(pair<hal_index_t, bool>(mPos, 
                                                         refDelta != mDelta));
    }
  }
}

void  MappedSegmentColCompareTest::compareArrays()
{
  CuAssertTrue(_testCase, _colArray.size() == _blockArray.size());

  double win = 0;
  double tot = 0;
  
  for (size_t i = 0; i < _colArray.size(); ++i)
  {
    map<hal_index_t, bool>& colEntry = _colArray[i];
    map<hal_index_t, bool>& blockEntry = _blockArray[i];
    
    if ( colEntry.size() == blockEntry.size())
       ++win;
    ++tot;
    continue;
    CuAssertTrue(_testCase, colEntry.size() == blockEntry.size());
    for (map<hal_index_t, bool>::iterator j = colEntry.begin(); 
         j != colEntry.end(); ++j)
    {
      map<hal_index_t, bool>::iterator k = blockEntry.find(j->first);
      CuAssertTrue(_testCase, k != blockEntry.end());
      CuAssertTrue(_testCase, k->second == j->second);
    }
  }
  cout << win / tot << endl;
}

void MappedSegmentColCompareTest1::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        0.75, 
                        0.1,
                        5,
                        10,
                        1000,
                        5,
                        10);
}

void MappedSegmentColCompareTest2::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        1.25, 
                        0.7,
                        10,
                        2,
                        50,
                        100,
                        5000);
}

void MappedSegmentColCompareTest3::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        1.25, 
                        0.7,
                        20,
                        2,
                        50,
                        1000,
                        7000);
}

void halMappedSegmentMapUpTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapUpTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapDownTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapDownTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapAcrossTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapAcrossTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapDupeTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapDupeTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentColCompareTest1(CuTest *testCase)
{
  try 
  {
    MappedSegmentColCompareTest1 tester;
    tester.check(testCase);
  }
  catch (double) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentColCompareTest2(CuTest *testCase)
{
  try 
  {
    MappedSegmentColCompareTest2 tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentColCompareTest3(CuTest *testCase)
{
  try 
  {
    MappedSegmentColCompareTest3 tester;
    tester.check(testCase);
  }
  catch (double) 
  {
    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halMappedSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halMappedSegmentMapUpTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDownTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapAcrossTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDupeTest);
  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest1);
  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest2);
  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest3);
  return suite;
}

