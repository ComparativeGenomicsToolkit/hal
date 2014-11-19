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
  // child1 and it is reversed and nonreversed to child2
  Genome* parent = alignment->addRootGenome("parent");
  Genome* child1 = alignment->addLeafGenome("child1", "parent", 1);
  Genome* child2 = alignment->addLeafGenome("child2", "parent", 1);
  // add a bunch of grandchildren with no rearrangemnts to test
  // simple parsing
  Genome* g1 = alignment->addLeafGenome("g1", "child2", 1);
  Genome* g2 = alignment->addLeafGenome("g2", "g1", 1);
  Genome* g3 = alignment->addLeafGenome("g3", "g2", 1);
  Genome* g4 = alignment->addLeafGenome("g4", "g3", 1);
  Genome* g5 = alignment->addLeafGenome("g5", "g4", 1);
  // add some with random inversions
  Genome* gi1 = alignment->addLeafGenome("gi1", "child1", 1);
  Genome* gi2 = alignment->addLeafGenome("gi2", "gi1", 1);
  Genome* gi3 = alignment->addLeafGenome("gi3", "gi2", 1);
  Genome* gi4 = alignment->addLeafGenome("gi4", "gi3", 1);
  Genome* gi5 = alignment->addLeafGenome("gi5", "gi4", 1);
  Genome* gs[] = {g1, g2, g3, g4, g5};
  Genome* gis[] = {gi1, gi2, gi3, gi4, gi5};
  seqVec[0] = Sequence::Info("Sequence", 12, 0, 1);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 1, 6);
  child1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 1, 6);
  child2->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 6, 4);
  g1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 4, 3);
  g2->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 3, 2);
  g3->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 2, 12);
  g4->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 12, 0);
  g5->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 6, 4);
  gi1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 4, 3);
  gi2->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 3, 2);
  gi3->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 2, 12);
  gi4->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 12, 12, 0);
  gi5->setDimensions(seqVec);


  parent->setString("CCCTACTTGTGC");
  child1->setString("CCCTACTTGTGC");
  child2->setString("CCCTACTTGTGC");
  for (size_t i = 0; i < 5; ++i)
  {
    gs[i]->setString("TCCTACTTGTGC");
    gis[i]->setString("TCCTACTTGTGC");
  }

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 12);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 12, 0, true, 0);
  ts.applyTo(ti);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 12, 0, false, 0);
  ts.applyTo(ti);
  
  for (size_t i = 0; i < 6; ++i)
  {
    bi = child2->getBottomSegmentIterator(i);
    bs.set(i * 2, 2, 0);
    bs._children.clear();
    bs._children.push_back(pair<hal_size_t, bool>(i, false));
    bs.applyTo(bi);

    ti = g1->getTopSegmentIterator(i);
    ts.set(i * 2, 2, i, false);
    ts.applyTo(ti);
  }

  for (size_t i = 0; i < 6; ++i)
  {
    bi = child1->getBottomSegmentIterator(i);
    bs.set(i * 2, 2, 0);
    bs._children.clear();
    bs._children.push_back(pair<hal_size_t, bool>(i, false));
    bs.applyTo(bi);

    ti = gi1->getTopSegmentIterator(i);
    ts.set(i * 2, 2, i, false);
    ts.applyTo(ti);
  }

  for (size_t i = 0; i < 5; ++i)
  {
    const Genome* g = gs[i];
    const Genome* parent = g->getParent();
    const Genome* child = i == 4 ? NULL : g->getChild(0);
    hal_size_t segLen = g->getSequenceLength() / g->getNumTopSegments();
    hal_size_t psegLen = parent->getSequenceLength() / 
       parent->getNumTopSegments();
    hal_size_t csegLen = 0;
    if (child)
    {
      csegLen =  child->getSequenceLength() / child->getNumTopSegments();
    }
    
    for (size_t j = 0; j < g->getNumTopSegments(); ++j)
    {
      bool inv = false;
      bi = parent->getBottomSegmentIterator(j);
      bs.set(j * segLen, segLen, (j * segLen) / psegLen);
      bs._children.clear();
      bs._children.push_back(pair<hal_size_t, bool>(j, inv));
      bs.applyTo(bi);

      hal_index_t bparse = NULL_INDEX;
      if (child != NULL)
      {
        bparse = (j * segLen) / csegLen;
      }
      ti = g->getTopSegmentIterator(j);
      ts.set(j * segLen, segLen, j, inv, bparse);
      ts.applyTo(ti);      
    }
  }
  
  for (size_t i = 0; i < 5; ++i)
  {
    const Genome* g = gis[i];
    const Genome* parent = g->getParent();
    const Genome* child = i == 4 ? NULL : g->getChild(0);
    hal_size_t segLen = g->getSequenceLength() / g->getNumTopSegments();
    hal_size_t psegLen = parent->getSequenceLength() / 
       parent->getNumTopSegments();
    hal_size_t csegLen = 0;
    if (child)
    {
      csegLen =  child->getSequenceLength() / child->getNumTopSegments();
    }
    
    for (size_t j = 0; j < g->getNumTopSegments(); ++j)
    {
      bool inv = rand() % 4 == 0;
      bi = parent->getBottomSegmentIterator(j);
      bs.set(j * segLen, segLen, (j * segLen) / psegLen);
      bs._children.clear();
      bs._children.push_back(pair<hal_size_t, bool>(j, inv));
      bs.applyTo(bi);

      hal_index_t bparse = NULL_INDEX;
      if (child != NULL)
      {
        bparse = (j * segLen) / csegLen;
      }
      ti = g->getTopSegmentIterator(j);
      ts.set(j * segLen, segLen, j, inv, bparse);
      ts.applyTo(ti);      
    }
  }

}

void MappedSegmentMapUpTest::testTopSegment(AlignmentConstPtr alignment,
                                            TopSegmentIteratorConstPtr top,
                                            const string& ancName)
{
  const Genome* parent = alignment->openGenome(ancName);
  set<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, parent, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = *results.begin();
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  // extra hop for when top is in grand child
  if (bottom->getGenome() != parent)
  {
    TopSegmentIteratorConstPtr temp = 
       bottom->getGenome()->getTopSegmentIterator();
    temp->toParseUp(bottom);
    bottom->toParent(temp);
  }
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
  validateAlignment(alignment);
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");
  TopSegmentIteratorConstPtr top = child2->getTopSegmentIterator();
  testTopSegment(alignment, top, "parent");
  top->slice(1,2);
  testTopSegment(alignment, top, "parent");
  top->toReverse();
  testTopSegment(alignment, top, "parent");
  top = child1->getTopSegmentIterator();
  testTopSegment(alignment, top, "parent");
  top->slice(1,2);
  testTopSegment(alignment, top, "parent");
  top->toReverse();
  testTopSegment(alignment, top, "parent");
  
  const Genome* g1 = alignment->openGenome("g1");
  for (hal_size_t i = 0; i < g1->getNumTopSegments(); ++i)
  {
    top = g1->getTopSegmentIterator(i);
    testTopSegment(alignment, top, "parent");
    top->slice(1,0);
    testTopSegment(alignment, top, "parent");
    top->toReverse();
    testTopSegment(alignment, top, "parent");
    top->slice(0,1);
    testTopSegment(alignment, top, "parent");
    top->toReverse();
    testTopSegment(alignment, top, "parent");
  }
}

void MappedSegmentParseTest::testTopSegment(AlignmentConstPtr alignment,
                                            TopSegmentIteratorConstPtr top)
{
  const Genome* parent = alignment->openGenome("parent");
  set<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, parent, NULL, false);

  vector<bool> covered(top->getLength(), false);
  
  CuAssertTrue(_testCase, results.size() >= 1);

  set<MappedSegmentConstPtr>::iterator i = results.begin();
  for (; i != results.end(); ++i)
  {
    MappedSegmentConstPtr mseg = *i;
    CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
    CuAssertTrue(_testCase, mseg->getGenome() == parent);
    for (hal_index_t j = mseg->getStartPosition(); j <= mseg->getEndPosition(); 
         ++j)
    {
      CuAssertTrue(_testCase, covered[j] == false);
      covered[j] = true;
    }
    CuAssertTrue(_testCase, mseg->getStartPosition() == 
                 mseg->getSource()->getStartPosition());
    CuAssertTrue(_testCase, mseg->getEndPosition() == 
                 mseg->getSource()->getEndPosition());

    set<MappedSegmentConstPtr> tResults;
    mseg->getMappedSegments(tResults, top->getGenome(), NULL, false);
    CuAssertTrue(_testCase, tResults.size() == 1);
    MappedSegmentConstPtr tmseg = *tResults.begin();
    CuAssertTrue(_testCase, tmseg->getGenome() == top->getGenome());
    CuAssertTrue(_testCase, tmseg->getSource()->getGenome() == 
                 mseg->getGenome());
    CuAssertTrue(_testCase, tmseg->getStartPosition() == 
                 mseg->getStartPosition());
    CuAssertTrue(_testCase, tmseg->getEndPosition() == 
                 mseg->getEndPosition());
    CuAssertTrue(_testCase, tmseg->getSource()->getStartPosition() == 
                 mseg->getStartPosition());
    CuAssertTrue(_testCase, tmseg->getSource()->getEndPosition() == 
                 mseg->getEndPosition());
  }
}

void MappedSegmentParseTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* g1 = alignment->openGenome("g1");
  const Genome* g2 = alignment->openGenome("g2");
  const Genome* g3 = alignment->openGenome("g3");
  const Genome* g4 = alignment->openGenome("g4");
  const Genome* g5 = alignment->openGenome("g5");
  const Genome* gs[] = {g1, g2, g3, g4, g5};
  
  for (size_t i = 0; i < 5; ++i)
  {
    const Genome* g = gs[i];
    for (size_t j = 0; j < g->getNumTopSegments(); ++j)
    {
      TopSegmentIteratorConstPtr top = g->getTopSegmentIterator(j);
      testTopSegment(alignment, top);
    }
  }
}

void MappedSegmentMapDownTest::testBottomSegment(
  AlignmentConstPtr alignment,
  BottomSegmentIteratorConstPtr bottom,
  hal_size_t childIndex)
{
  const Genome* child = bottom->getGenome()->getChild(childIndex);
  set<MappedSegmentConstPtr> results;
  bottom->getMappedSegments(results, child, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = *results.begin();
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
  validateAlignment(alignment);
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
  set<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, other, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = *results.begin();
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
  validateAlignment(alignment);
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
  bs.set(0, 3);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 3, 0, true, NULL_INDEX, 1);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, true, NULL_INDEX, 2);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, 0, true, NULL_INDEX, 0);
  ts.applyTo(ti);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 3, 0, false);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, NULL_INDEX, true);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, NULL_INDEX, false);
  ts.applyTo(ti);
}

void MappedSegmentMapDupeTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");

  TopSegmentIteratorConstPtr top = child1->getTopSegmentIterator();
  set<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, child2, NULL, true);
//  CuAssertTrue(_testCase, results.size() == 3);
  
  MappedSegmentConstPtr mseg = *results.begin();
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
  set<MappedSegmentConstPtr>::iterator i = results.begin();
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

void MappedSegmentMapExtraParalogsTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);

  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  // Set up a case where all the segments of grandChild1 coalesce with
  // the first segment of grandChild2, but only if using the root as
  // the coalescence limit. Otherwise only the first segments map to
  // each other.
  Genome* root = alignment->addRootGenome("root");
  Genome* parent = alignment->addLeafGenome("parent", "root", 1);
  Genome* grandChild1 = alignment->addLeafGenome("grandChild1", "parent", 1);
  Genome* grandChild2 = alignment->addLeafGenome("grandChild2", "parent", 1);
  seqVec[0] = Sequence::Info("Sequence", 3, 0, 1);
  root->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 3);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  grandChild1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  grandChild2->setDimensions(seqVec);

  root->setString("CCC");
  parent->setString("CCCTACGTG");
  grandChild1->setString("CCCTACGTG");
  grandChild2->setString("CCCTACGTG");

  bi = root->getBottomSegmentIterator();
  bs.set(0, 3);
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);

  ti = parent->getTopSegmentIterator();
  ts.set(0, 3, 0, false, NULL_INDEX, 1);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, false, NULL_INDEX, 2);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, 0, false, NULL_INDEX, 0);
  ts.applyTo(ti);

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 3);
  bs._children.clear();
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
  bi->toRight();
  bs.set(3, 3);
  bs._children.clear();
  bs._children.push_back(pair<hal_size_t, bool>(1, true));
  bs._children.push_back(pair<hal_size_t, bool>(NULL_INDEX, true));
  bs.applyTo(bi);
  bi->toRight();
  bs.set(6, 3);
  bs._children.clear();
  bs._children.push_back(pair<hal_size_t, bool>(2, true));
  bs._children.push_back(pair<hal_size_t, bool>(NULL_INDEX, false));
  bs.applyTo(bi);

  ti = grandChild1->getTopSegmentIterator();
  ts.set(0, 3, 0, true);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 1, true);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, 2, true);
  ts.applyTo(ti);

  ti = grandChild2->getTopSegmentIterator();
  ts.set(0, 3, 0, false);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, NULL_INDEX, true);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(6, 3, NULL_INDEX, false);
  ts.applyTo(ti);

  parent->fixParseInfo();
}

void MappedSegmentMapExtraParalogsTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome *grandChild1 = alignment->openGenome("grandChild1");
  const Genome *grandChild2 = alignment->openGenome("grandChild2");
  const Genome *root = alignment->openGenome("root");

  TopSegmentIteratorConstPtr top = grandChild2->getTopSegmentIterator();
  set<MappedSegmentConstPtr> results;

  // First, check that by default we will only get the homologies in
  // or before the MRCA. (in this case, just seg 0 of grandChild1).
  top->getMappedSegments(results, grandChild1, NULL, true);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = *results.begin();
  // Source information should be preserved
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());

  // Check target information is correct
  CuAssertTrue(_testCase,
               mseg->getGenome() == grandChild1);
  CuAssertTrue(_testCase,
               mseg->getStartPosition() == 2);
  CuAssertTrue(_testCase,
               mseg->getLength() == 3);
  CuAssertTrue(_testCase,
               mseg->getReversed() == true);

  // Check that by using the grandparent as the coalescence limit we
  // will get all the paralogs.
  top->getMappedSegments(results, grandChild1, NULL, true, 0, root);
  CuAssertTrue(_testCase, results.size() == 3);
  set<MappedSegmentConstPtr>::iterator i = results.begin();
  bool found[3] = {false, false, false};
  for (; i != results.end(); ++i)
  {
      // Source information should be preserved
    CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
    CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
                 top->getStartPosition());
    CuAssertTrue(_testCase, 
                 mseg->getSource()->getLength() == top->getLength());
    CuAssertTrue(_testCase, 
                 mseg->getSource()->getReversed() == top->getReversed());
    
    // Check target information is correct
    CuAssertTrue(_testCase,
                 mseg->getGenome() == grandChild1);
    CuAssertTrue(_testCase,
                 mseg->getStartPosition() == 2
                 || mseg->getStartPosition() == 5
                 || mseg->getStartPosition() == 8);
    CuAssertTrue(_testCase,
                 mseg->getLength() == 3);
    CuAssertTrue(_testCase,
                 mseg->getReversed() == true);
    found[mseg->getArrayIndex()] = true;
  }
}

void  MappedSegmentColCompareTest::checkCallBack(AlignmentConstPtr alignment)
{
  if (alignment->getNumGenomes() == 0)
  {
    return;
  }

  validateAlignment(alignment);
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
        _ref = srcGenome;
        _tgt = tgtGenome;
        createColArray();
        createBlockArray();
        compareArrays();
      }
    }
  }
}

void  MappedSegmentColCompareTest::createColArray()
{
  hal_size_t N = _ref->getSequenceLength();
  _colArray.clear();
  _colArray.resize(N);
  set<const Genome*> tgtSet;
  tgtSet.insert(_tgt);
  ColumnIteratorConstPtr colIt = _ref->getColumnIterator(&tgtSet, 0, 0, 
                                                         NULL_INDEX, false,
                                                         false, false, true);
  while (true)
  {
    const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
    ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
    vector<pair<hal_index_t, bool> > insertList;
    // Pass 1 find all homologies in target
    for (; colMapIt != colMap->end(); colMapIt++)
    {
      if (colMapIt->first->getGenome() == _tgt)
      {
        ColumnIterator::DNASet* dnaSet = colMapIt->second;
        ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
        for (; dnaIt != dnaSet->end(); ++dnaIt)
        {
          DNAIteratorConstPtr dna = *dnaIt;
          insertList.push_back(
            pair<hal_index_t, bool>(dna->getArrayIndex(), dna->getReversed()));
        }
      }
      else
      {
        CuAssertTrue(_testCase, colMapIt->first->getGenome() == _ref);
      }
    }

    // Pass 2 update each reference position with all homologies found
    // in Pass 1
    for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++)
    {
      if (colMapIt->first->getGenome() == _ref)
      {
        ColumnIterator::DNASet* dnaSet = colMapIt->second;
        ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
        for (; dnaIt != dnaSet->end(); ++dnaIt)
        {
          DNAIteratorConstPtr dna = *dnaIt;
          for (size_t insIdx = 0; insIdx < insertList.size(); ++insIdx)
          {
            _colArray[dna->getArrayIndex()].insert(insertList[insIdx]);
          }
        }
      }
    }   
    if (colIt->lastColumn())
    {
      break;
    }
    colIt->toRight();
  }
}

void  MappedSegmentColCompareTest::createBlockArray()
{
  hal_size_t N = _ref->getSequenceLength();
  _blockArray.clear();
  _blockArray.resize(N);
  SegmentIteratorConstPtr refSeg;
  hal_index_t numSegs;
  if (_ref->getNumTopSegments() > 0)
  { 
    refSeg = _ref->getTopSegmentIterator(0);
    numSegs = _ref->getNumTopSegments();
  }
  else
  {
    refSeg = _ref->getBottomSegmentIterator(0);
    numSegs = _ref->getNumBottomSegments();
  }

  set<MappedSegmentConstPtr> results;  
  for (; refSeg->getArrayIndex() < numSegs; refSeg->toRight())
  {
    refSeg->getMappedSegments(results, _tgt);
  }
  
  for (set<MappedSegmentConstPtr>::iterator i = results.begin();
       i != results.end(); ++i)
  {
    MappedSegmentConstPtr mseg = *i;
    CuAssertTrue(_testCase, mseg->getLength() == mseg->getSource()->getLength());
    SlicedSegmentConstPtr refSeg = mseg->getSource();
    hal_index_t refDelta = refSeg->getReversed() ? -1 : 1;
    CuAssertTrue(_testCase, refDelta == 1);
    hal_index_t mDelta = mseg->getReversed() ? -1 : 1;
    if (refSeg->getGenome() != _ref)
       cout << "ref " << refSeg->getGenome()->getName() << " != " 
            << _ref->getName() << " note taget is "
            << mseg->getGenome()->getName() << endl;
    CuAssertTrue(_testCase, refSeg->getGenome() == _ref);

    for (hal_index_t offset = 0; offset < (hal_index_t)mseg->getLength(); 
         ++offset)
    {
      hal_index_t refPos = refSeg->getStartPosition() + offset * refDelta;
      hal_index_t mPos = mseg->getStartPosition() + offset * mDelta;
      _blockArray[refPos].insert(pair<hal_index_t, bool>(mPos, 
                                                         refDelta != mDelta));
    }
  }
}

void MappedSegmentColCompareTest::compareArrays()
{
  CuAssertTrue(_testCase, _colArray.size() == _blockArray.size());

  for (size_t i = 0; i < _colArray.size(); ++i)
  {
    set<pair<hal_index_t, bool> >& colEntry = _colArray[i];
    set<pair<hal_index_t, bool> >& blockEntry = _blockArray[i];
    
    if (colEntry.size() != blockEntry.size())
    {
      cout << "col " << i << " cs=" << colEntry.size() << " bs="
           << blockEntry.size() << endl;
      //continue;
    }
    CuAssertTrue(_testCase, colEntry.size() == blockEntry.size());
    for (set<pair<hal_index_t, bool> >::iterator j = colEntry.begin(); 
         j != colEntry.end(); ++j)
    {
      set<pair<hal_index_t, bool> >::iterator k = blockEntry.find(*j);

      // because of global dupe-skipping built into the column iterator
      // and not presently in (or wanted?) in blockmapper the orientation of 
      // paralogies can be flipped wrt the reference.  This means that
      // any homolgoy can be flipped if it went through a dupe cycle so
      // we can't compare inversion flags
      set<pair<hal_index_t, bool> >::iterator kk = blockEntry.end();
      kk = blockEntry.find(pair<hal_index_t, bool>(j->first, !j->second));
      if (k == blockEntry.end() && kk == blockEntry.end())
      {
        cout << j->first << ", " << j->second << " not found for ref " << i 
             << endl;
      }
      CuAssertTrue(_testCase, k != blockEntry.end() || kk != blockEntry.end());
    }
  }
}

void MappedSegmentColCompareTestCheck1::createCallBack(AlignmentPtr alignment)
{
  MappedSegmentMapUpTest::createCallBack(alignment);
}

void 
MappedSegmentColCompareTestCheck1::checkCallBack(AlignmentConstPtr alignment)
{
  MappedSegmentColCompareTest::checkCallBack(alignment);
}

void MappedSegmentColCompareTestCheck2::createCallBack(AlignmentPtr alignment)
{
  MappedSegmentMapDupeTest::createCallBack(alignment);
}

void 
MappedSegmentColCompareTestCheck2::checkCallBack(AlignmentConstPtr alignment)
{
  MappedSegmentColCompareTest::checkCallBack(alignment);
}

void MappedSegmentColCompareTest1::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        2, // meanDegree
                        0.1, // maxBranchLength
                        6, // maxGenomes
                        10, // minSegmentLength
                        1000, // maxSegmentLength
                        5, // minSegments
                        10, // maxSegments
                        1101); // seed
}

void MappedSegmentColCompareTest2::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        1.25, 
                        0.7,
                        8,
                        2,
                        50,
                        10,
                        500,
                        23);
}

void MappedSegmentColCompareTest3::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        1.5, 
                        0.7,
                        12,
                        1,
                        100,
                        50,
                        1000);
}

void halMappedSegmentMapUpTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapUpTest tester;
    tester.check(testCase);
  }
  catch (double) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentParseTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentParseTest tester;
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

void halMappedSegmentMapExtraParalogsTest(CuTest *testCase)
{
/*  try 
  {*/
    MappedSegmentMapExtraParalogsTest tester;
    tester.check(testCase);
/*  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } */
}

void haMappedSegmentColCompareTestCheck1(CuTest *testCase)
{
  try 
  {
    MappedSegmentColCompareTestCheck1 tester;
    tester.check(testCase);
  }
  catch (double) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void haMappedSegmentColCompareTestCheck2(CuTest *testCase)
{
  try 
  {
    MappedSegmentColCompareTestCheck2 tester;
    tester.check(testCase);
  }
  catch (double) 
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
  SUITE_ADD_TEST(suite, halMappedSegmentMapExtraParalogsTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapUpTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDownTest);
  SUITE_ADD_TEST(suite, halMappedSegmentParseTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapAcrossTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDupeTest); 
  SUITE_ADD_TEST(suite, haMappedSegmentColCompareTestCheck1);
  SUITE_ADD_TEST(suite, haMappedSegmentColCompareTestCheck2);
  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest1);
//  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest2);
//  SUITE_ADD_TEST(suite, halMappedSegmentColCompareTest3);
  return suite;
}

