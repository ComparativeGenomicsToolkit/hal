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
#include "halColumnIteratorTest.h"
#include "halRandomData.h"
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"
#include "hal.h"


using namespace std;
using namespace hal;

void ColumnIteratorBaseTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 10, 1e-10, 3, 77, 77, 10, 10);
}

void ColumnIteratorBaseTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome(alignment->getRootName());
  const Sequence* sequence = genome->getSequenceBySite(0);
  assert(genome != NULL);

  // Iterate over the genome, ensuring that base i aligns to single 
  // other base
  ColumnIteratorConstPtr colIterator = sequence->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    CuAssertTrue(_testCase, colMap->size() == 3);
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      ColumnIterator::DNASet::const_iterator dnaIt = i->second->begin();
      for (size_t j = 0; j < i->second->size(); ++j)
      {
        CuAssertTrue(_testCase, i->second->size() == 1);
        CuAssertTrue(_testCase, (*dnaIt)->getArrayIndex() == 
                     (hal_index_t)columnNumber);
        dnaIt++;
      }
    }

    colIterator->toRight();
  }
}

void ColumnIteratorDepthTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("grandpa");
  alignment->addLeafGenome("dad", "grandpa", branchLength);
  alignment->addLeafGenome("son1", "dad", branchLength);
  alignment->addLeafGenome("son2", "dad", branchLength);
  
  vector<Sequence::Info> dims(1);
  hal_size_t numSegments = 10;
  hal_size_t segLength = 10;
  hal_size_t seqLength = numSegments * segLength;

  Genome* son1 = alignment->openGenome("son1");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son1->setDimensions(dims);
  son1->setString(string(seqLength, 'A'));
  
  Genome* son2 = alignment->openGenome("son2");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son2->setDimensions(dims);
  son2->setString(string(seqLength, 'C'));

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, numSegments);
  dad->setDimensions(dims);
  dad->setString(string(seqLength, 'G'));

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("seq", seqLength, 0, numSegments);
  grandpa->setDimensions(dims);
  grandpa->setString(string(seqLength, 'T'));

  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = son2->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = dad->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i, false, i);
    ts.applyTo(ti);

    bi = dad->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildReversed(0, false);
    bi->getBottomSegment()->setChildIndex(1, i);
    bi->getBottomSegment()->setChildReversed(1, false);;
    bs.set(i * segLength, segLength, i);
    bs.applyTo(bi);

    bi = grandpa->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildReversed(0, false);
    bs.set(i * segLength, segLength);
    bs.applyTo(bi);
  }  
}

void ColumnIteratorDepthTest::checkGenome(const Genome* genome)
{
  assert(genome != NULL);
  const Sequence* sequence = genome->getSequenceBySite(0);
  ColumnIteratorConstPtr colIterator = sequence->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    CuAssertTrue(_testCase, colMap->size() == 4);
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      CuAssertTrue(_testCase, i->second->size() == 1);
      DNAIteratorConstPtr dnaIt = *i->second->begin();
      
/*        cout << "column=" << columnNumber 
           << " genome=" << dnaIt->getGenome()->getName()
           << " index=" << dnaIt->getArrayIndex() << endl;
*/
      CuAssertTrue(_testCase, dnaIt->getArrayIndex() == 
                   (hal_index_t)columnNumber);
    }
    colIterator->toRight();
  }
}

void ColumnIteratorDepthTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome("grandpa");
  checkGenome(genome);
  genome = alignment->openGenome("dad");
  checkGenome(genome);
  genome = alignment->openGenome("son1");
  checkGenome(genome);
  genome = alignment->openGenome("son2");
  checkGenome(genome);
}

void ColumnIteratorDupTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("dad");
  alignment->addLeafGenome("son1", "dad", branchLength);
  alignment->addLeafGenome("son2", "dad", branchLength);
  
  vector<Sequence::Info> dims(1);
  hal_size_t numSegments = 10;
  hal_size_t segLength = 10;
  hal_size_t seqLength = numSegments * segLength;

  Genome* son1 = alignment->openGenome("son1");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son1->setDimensions(dims);
  son1->setString(string(seqLength, 't'));
  
  Genome* son2 = alignment->openGenome("son2");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son2->setDimensions(dims);
  son2->setString(string(seqLength, 'c'));

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("seq", seqLength, 0, numSegments);
  dad->setDimensions(dims);
  dad->setString(string(seqLength, 'G'));
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = son2->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    bi = dad->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildReversed(0, false);
    bi->getBottomSegment()->setChildIndex(1, i);
    bi->getBottomSegment()->setChildReversed(1, false);
    bs.set(i * segLength, segLength);
    bs.applyTo(bi);
  }  
  
  // son1 is just one big duplication stemming 
  // from 0th segment in dad
  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ti->getTopSegment()->setParentIndex(0);
    if (i < numSegments - 1)
    {
      ti->getTopSegment()->setNextParalogyIndex(i+1);
    }
    else
    {
      ti->getTopSegment()->setNextParalogyIndex(0);
    }
    if (i > 0)
    {
      bi = dad->getBottomSegmentIterator(i);
      bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
    }
  }
  
  // son2 has a duplication between 4 and 8,
  // which both derive from 4 in dad
  ti = son2->getTopSegmentIterator(4);
  ti->getTopSegment()->setNextParalogyIndex(8);
  ti = son2->getTopSegmentIterator(8);
  ti->getTopSegment()->setNextParalogyIndex(4);
  ti->getTopSegment()->setParentIndex(4);
  bi = dad->getBottomSegmentIterator(8);
  bi->getBottomSegment()->setChildIndex(1, NULL_INDEX);
}

void ColumnIteratorDupTest::checkGenome(const Genome* genome)
{
  assert(genome != NULL);
  const Sequence* sequence = genome->getSequenceBySite(0);
  ColumnIteratorConstPtr colIterator = sequence->getColumnIterator();
  size_t colNumber = 0;
  for (; colNumber < genome->getSequenceLength(); colIterator->toRight(), 
          ++colNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    // check that all three genomes are in the map
    CuAssertTrue(_testCase, colMap->size() == 3);

    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      DNAIteratorConstPtr dnaIt = *i->second->begin();
      // the first segment (of any genome) should be aligned to 
      // every segment in son1
      if (i->first->getGenome()->getName() == "son1" && 
          genome->getName() != "son1")
      {
        if (colNumber < i->first->getTopSegmentIterator()->getLength())
        {
          CuAssertTrue(_testCase, i->second->size() ==
                       i->first->getTopSegmentIterator()->getLength());
        }
        else
        {
          CuAssertTrue(_testCase, i->second->size() == 0);
        }
      }
      // check the paralogy on son2 
      else if (i->first->getGenome()->getName() == "son2" && genome->getName() == "dad")
      {
        if (colNumber >= 40 && colNumber < 50)
        {
          CuAssertTrue(_testCase, i->second->size() == 2);
        }
        else if (colNumber >= 80 && colNumber < 90)
        {
          CuAssertTrue(_testCase, i->second->size() == 0);
        }
        else
        {
          CuAssertTrue(_testCase, i->second->size() == 1);
        }
      }
    }
  }
}

void ColumnIteratorDupTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome("dad");
  checkGenome(genome);
  genome = alignment->openGenome("son1");
  checkGenome(genome);
  genome = alignment->openGenome("son2");
  checkGenome(genome);
}

void ColumnIteratorInvTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("grandpa");
  alignment->addLeafGenome("dad", "grandpa", branchLength);
  alignment->addLeafGenome("son1", "dad", branchLength);
  
  vector<Sequence::Info> dims(1);
  hal_size_t numSegments = 10;
  hal_size_t segLength = 10;
  hal_size_t seqLength = numSegments * segLength;

  Genome* son1 = alignment->openGenome("son1");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son1->setDimensions(dims);

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, numSegments);
  dad->setDimensions(dims);

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("seq", seqLength, 0, numSegments);
  grandpa->setDimensions(dims);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);
    ti = dad->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i, false, i);
    ts.applyTo(ti);

    bi = dad->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildReversed(0, false);
    bs.set(i * segLength, segLength, i);
    bs.applyTo(bi);

    bi = grandpa->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildReversed(0, false);
    bs.set(i * segLength, segLength);
    bs.applyTo(bi);
  }  

  DNAIteratorPtr gDnaIt = grandpa->getDNAIterator();
  DNAIteratorPtr dDnaIt = dad->getDNAIterator();
  DNAIteratorPtr sDnaIt = son1->getDNAIterator();
  
  for (hal_size_t i = 0; i < grandpa->getSequenceLength(); ++i)
  {
    switch (i % 4) {
    case 0: 
      gDnaIt->setChar('A'); 
      dDnaIt->setChar('C'); 
      sDnaIt->setChar('T');
      break;
    case 1: 
      gDnaIt->setChar('G'); 
      dDnaIt->setChar('A'); 
      sDnaIt->setChar('C');
      break;
    case 2: 
      gDnaIt->setChar('T'); 
      dDnaIt->setChar('G'); 
      sDnaIt->setChar('A');
      break;
    case 3: 
      gDnaIt->setChar('C'); 
      dDnaIt->setChar('T'); 
      sDnaIt->setChar('G');
      break;
    }
    gDnaIt->toRight();
    dDnaIt->toRight();
    sDnaIt->toRight();
  }

  // child-dad edge has inversion in 0th segment
  bi = dad->getBottomSegmentIterator();
  bi->getBottomSegment()->setChildReversed(0, true);
  ti = son1->getTopSegmentIterator();
  ti->getTopSegment()->setParentReversed(true);

  // child-dad edge and dad-grandpa edges have inversion
  // in 1st segment
  bi = dad->getBottomSegmentIterator(1);
  bi->getBottomSegment()->setChildReversed(0, true);
  ti = son1->getTopSegmentIterator(1);
  ti->getTopSegment()->setParentReversed(true);
  bi = grandpa->getBottomSegmentIterator(1);
  bi->getBottomSegment()->setChildReversed(0, true);
  ti = dad->getTopSegmentIterator(1);
  ti->getTopSegment()->setParentReversed(true);
}

void ColumnIteratorInvTest::checkGenome(const Genome* genome)
{
  assert(genome != NULL);
  const Sequence* sequence = genome->getSequenceBySite(0);
  ColumnIteratorConstPtr colIterator = sequence->getColumnIterator();
  size_t colNumber = 0;
  for (; colNumber < genome->getSequenceLength(); colIterator->toRight(), 
          ++colNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();

    // check that all three sequences are in the map
    CuAssertTrue(_testCase, colMap->size() == 3);

    // iterate the sequences in the map
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      CuAssertTrue(_testCase, i->second->size() == 1);
      DNAIteratorConstPtr dnaIt = *i->second->begin();
      
      if (i->first->getGenome()->getName() == "son1")
      {
        const Genome* dad = i->first->getGenome()->getParent();
        const Genome* grandpa = dad->getParent();

        CuAssertTrue(_testCase, dad->getName() == "dad");
        CuAssertTrue(_testCase, grandpa->getName() == "grandpa");
        
        SequenceIteratorConstPtr dadSeq = dad->getSequenceIterator();
        SequenceIteratorConstPtr graSeq = grandpa->getSequenceIterator();
        
        CuAssertTrue(_testCase, 
                     dadSeq->getSequence()->getGenome() == dad);
        CuAssertTrue(_testCase, 
                     graSeq->getSequence()->getGenome() == grandpa);
        
        CuAssertTrue(_testCase, colMap->find(dadSeq->getSequence()) != 
                     colMap->end());
        CuAssertTrue(_testCase, colMap->find(graSeq->getSequence()) != 
                     colMap->end());

        DNAIteratorConstPtr dadIt = 
           *colMap->find(dadSeq->getSequence())->second->begin();
        DNAIteratorConstPtr graIt = 
           *colMap->find(graSeq->getSequence())->second->begin();

        if (colNumber >= 0 && colNumber < 10)
        {
          CuAssertTrue(_testCase, dnaIt->getReversed() == false);
          CuAssertTrue(_testCase, 
                       dnaIt->getArrayIndex() == (hal_index_t)colNumber); 

          // inversion on bottom branch
          CuAssertTrue(_testCase, dadIt->getReversed() == true);
          CuAssertTrue(_testCase, 
                       dadIt->getArrayIndex() == 9 - (hal_index_t)colNumber); 

          CuAssertTrue(_testCase, graIt->getReversed() == true);
          CuAssertTrue(_testCase, 
                       graIt->getArrayIndex() == 9 - (hal_index_t)colNumber); 
        }
        else if (colNumber >= 10 && colNumber < 20)
        {
          CuAssertTrue(_testCase, dnaIt->getReversed() == false);
          CuAssertTrue(_testCase, 
                       dnaIt->getArrayIndex() == (hal_index_t)colNumber);
 
          // inversion on bottom branch
          CuAssertTrue(_testCase, dadIt->getReversed() == true);
          CuAssertTrue(_testCase, 
                       dadIt->getArrayIndex() == 29 - (hal_index_t)colNumber);
 
          // inversion on top branch
          CuAssertTrue(_testCase, graIt->getReversed() == false);
          CuAssertTrue(_testCase, 
                       graIt->getArrayIndex() == (hal_index_t)colNumber); 
        }
        else
        {
          CuAssertTrue(_testCase, 
                       dnaIt->getArrayIndex() == (hal_index_t)colNumber);
          CuAssertTrue(_testCase, 
                       dadIt->getArrayIndex() == (hal_index_t)colNumber);
          CuAssertTrue(_testCase, 
                       graIt->getArrayIndex() == (hal_index_t)colNumber);
        }         
      }
    }
  }
}

void ColumnIteratorInvTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome("son1");
  checkGenome(genome);
}

void ColumnIteratorGapTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("grandpa");
  alignment->addLeafGenome("dad", "grandpa", branchLength);

  vector<Sequence::Info> dims(1);

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("gseq", 12, 0, 3);
  grandpa->setDimensions(dims);

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("dseq", 8, 2, 0);
  dad->setDimensions(dims);

  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  // seg - seg - seg
  // seg - GAP - seg
  bi = grandpa->getBottomSegmentIterator(0);
  bs.set(0, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 0);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = grandpa->getBottomSegmentIterator(1);
  bs.set(4, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);

  bi = grandpa->getBottomSegmentIterator(2);
  bs.set(8, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 1);
  bi->getBottomSegment()->setChildReversed(0, false);

  ti = dad->getTopSegmentIterator(0);
  ts.set(0, 4, 0);
  ts.applyTo(ti);

  ti = dad->getTopSegmentIterator(1);
  ts.set(4, 4, 2);
  ts.applyTo(ti);

  grandpa->setString("ACGTAAAAGGGG");
  dad->setString("ACGTGGGG");
}

void ColumnIteratorGapTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* dad = alignment->openGenome("dad");
  const Sequence* dadSeq = dad->getSequence("dseq");
  const Genome* grandpa = alignment->openGenome("grandpa");
  const Sequence* grandpaSeq = grandpa->getSequence("gseq");

  ColumnIteratorConstPtr colIterator = dadSeq->getColumnIterator(NULL, 1000);
  for (hal_index_t i = 0; i < 12; ++i)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    if (i < 4 || i >= 8)
    {
      hal_index_t dadCol = i < 4 ? i : i - 4;
      ColumnIterator::DNASet* entry = colMap->find(dadSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);    
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == dadSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == dadCol);
    }

    ColumnIterator::DNASet* entry = colMap->find(grandpaSeq)->second;
    CuAssertTrue(_testCase, entry->size() == 1);    
    DNAIteratorConstPtr dna = entry->at(0);
    CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
    CuAssertTrue(_testCase, dna->getArrayIndex() == i);

    colIterator->toRight();
  }
}

void ColumnIteratorMultiGapTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("adam");
  alignment->addLeafGenome("grandpa", "adam", branchLength);
  alignment->addLeafGenome("dad", "grandpa", branchLength);

  vector<Sequence::Info> dims(1);

  Genome* adam = alignment->openGenome("adam");
  dims[0] = Sequence::Info("aseq", 16, 0, 4);
  adam->setDimensions(dims);

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("gseq", 12, 3, 3);
  grandpa->setDimensions(dims);

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("dseq", 8, 2, 0);
  dad->setDimensions(dims);

  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  // seg - seg - GAP - seg
  // seg - seg - seg
  // seg = GAP - seg
  bi = adam->getBottomSegmentIterator(0);
  bs.set(0, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 0);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = adam->getBottomSegmentIterator(1);
  bs.set(4, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 1);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = adam->getBottomSegmentIterator(2);
  bs.set(8, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  
  bi = adam->getBottomSegmentIterator(3);
  bs.set(12, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 2);
  bi->getBottomSegment()->setChildReversed(0, false);

  ti = grandpa->getTopSegmentIterator(0);
  ts.set(0, 4, 0, false, 0);
  ts.applyTo(ti);

  ti = grandpa->getTopSegmentIterator(1);
  ts.set(4, 4, 1, false, 1);
  ts.applyTo(ti);

  ti = grandpa->getTopSegmentIterator(2);
  ts.set(8, 4, 3, false, 2);
  ts.applyTo(ti);

  bi = grandpa->getBottomSegmentIterator(0);
  bs.set(0, 4, 0);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 0);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = grandpa->getBottomSegmentIterator(1);
  bs.set(4, 4, 1);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = grandpa->getBottomSegmentIterator(2);
  bs.set(8, 4, 2);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 1);
  bi->getBottomSegment()->setChildReversed(0, false);

  ti = dad->getTopSegmentIterator(0);
  ts.set(0, 4, 0);
  ts.applyTo(ti);

  ti = dad->getTopSegmentIterator(1);
  ts.set(4, 4, 2);
  ts.applyTo(ti);

  adam->setString("ACGTAAAATTTTGGGG");
  grandpa->setString("ACGTAAAAGGGG");
  dad->setString("ACGTGGGG");
}

void ColumnIteratorMultiGapTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);

  const Genome* dad = alignment->openGenome("dad");
  const Sequence* dadSeq = dad->getSequence("dseq");
  const Genome* grandpa = alignment->openGenome("grandpa");
  const Sequence* grandpaSeq = grandpa->getSequence("gseq");
  const Genome* adam = alignment->openGenome("adam");
  const Sequence* adamSeq = adam->getSequence("aseq");


  BottomSegmentIteratorConstPtr bi = grandpa->getBottomSegmentIterator(2);
  TopSegmentIteratorConstPtr ti = dad->getTopSegmentIterator(1);
  
  assert (bi->getBottomSegment()->getChildReversed(0) == 
          ti->getTopSegment()->getParentReversed());
  
  CuAssertTrue(_testCase,
               dad && dadSeq && grandpa && grandpaSeq && adam && adamSeq);

  ColumnIteratorConstPtr colIterator = dadSeq->getColumnIterator(NULL, 1000);

  for (hal_index_t i = 0; i < 16; ++i)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();

    if (i < 4)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(dadSeq) != colMap->end());
      entry = colMap->find(dadSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == dadSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
    }
    else if (i < 8)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
      
      CuAssertTrue(_testCase, colMap->find(dadSeq) == colMap->end() ||
                   colMap->find(dadSeq)->second->size() == 0);      
    }
    else if (i < 12)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);

      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
      
      CuAssertTrue(_testCase, colMap->find(grandpaSeq) == colMap->end() ||
                   colMap->find(grandpaSeq)->second->size() == 0);      

      CuAssertTrue(_testCase, colMap->find(dadSeq) == colMap->end() ||
                   colMap->find(dadSeq)->second->size() == 0);      
    }
    else
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i - 4);

      CuAssertTrue(_testCase, colMap->find(dadSeq) != colMap->end());      
      entry = colMap->find(dadSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == dadSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i - 8);
    }

    colIterator->toRight();
  }
}

void ColumnIteratorMultiGapInvTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("adam");
  alignment->addLeafGenome("grandpa", "adam", branchLength);
  alignment->addLeafGenome("dad", "grandpa", branchLength);

  vector<Sequence::Info> dims(1);

  Genome* adam = alignment->openGenome("adam");
  dims[0] = Sequence::Info("aseq", 16, 0, 4);
  adam->setDimensions(dims);

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("gseq", 12, 3, 3);
  grandpa->setDimensions(dims);

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("dseq", 8, 2, 0);
  dad->setDimensions(dims);

  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  // seg - seg - GAP - seg
  // seg - seg - seg
  // seg = GAP - seg
  bi = adam->getBottomSegmentIterator(0);
  bs.set(0, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 0);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = adam->getBottomSegmentIterator(1);
  bs.set(4, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 1);
  bi->getBottomSegment()->setChildReversed(0, true);

  bi = adam->getBottomSegmentIterator(2);
  bs.set(8, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  
  bi = adam->getBottomSegmentIterator(3);
  bs.set(12, 4);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 2);
  bi->getBottomSegment()->setChildReversed(0, false);

  ti = grandpa->getTopSegmentIterator(0);
  ts.set(0, 4, 0, false, 0);
  ts.applyTo(ti);

  ti = grandpa->getTopSegmentIterator(1);
  ts.set(4, 4, 1, true, 1);
  ts.applyTo(ti);

  ti = grandpa->getTopSegmentIterator(2);
  ts.set(8, 4, 3, false, 2);
  ts.applyTo(ti);

  bi = grandpa->getBottomSegmentIterator(0);
  bs.set(0, 4, 0);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 0);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = grandpa->getBottomSegmentIterator(1);
  bs.set(4, 4, 1);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  bi->getBottomSegment()->setChildReversed(0, false);

  bi = grandpa->getBottomSegmentIterator(2);
  bs.set(8, 4, 2);
  bs.applyTo(bi);
  bi->getBottomSegment()->setChildIndex(0, 1);
  bi->getBottomSegment()->setChildReversed(0, false);

  ti = dad->getTopSegmentIterator(0);
  ts.set(0, 4, 0);
  ts.applyTo(ti);

  ti = dad->getTopSegmentIterator(1);
  ts.set(4, 4, 2);
  ts.applyTo(ti);

  adam->setString("ACGTAAAATTTTGGGG");
  grandpa->setString("ACGTAAAAGGGG");
  dad->setString("ACGTGGGG");
}

void ColumnIteratorMultiGapInvTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* dad = alignment->openGenome("dad");
  const Sequence* dadSeq = dad->getSequence("dseq");
  const Genome* grandpa = alignment->openGenome("grandpa");
  const Sequence* grandpaSeq = grandpa->getSequence("gseq");
  const Genome* adam = alignment->openGenome("adam");
  const Sequence* adamSeq = adam->getSequence("aseq");

  CuAssertTrue(_testCase,
               dad && dadSeq && grandpa && grandpaSeq && adam && adamSeq);

  ColumnIteratorConstPtr colIterator = dadSeq->getColumnIterator(NULL, 1000);
  for (hal_index_t i = 0; i < 16; ++i)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();

    if (i < 4)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(dadSeq) != colMap->end());
      entry = colMap->find(dadSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == dadSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
    }
    else if (i < 8)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == 11 - i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
      
      CuAssertTrue(_testCase, colMap->find(dadSeq) == colMap->end() ||
                   colMap->find(dadSeq)->second->size() == 0);      
    }
    else if (i < 12)
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);

      CuAssertTrue(_testCase, dna->getArrayIndex() == i);
      
      CuAssertTrue(_testCase, colMap->find(grandpaSeq) == colMap->end() ||
                   colMap->find(grandpaSeq)->second->size() == 0);      

      CuAssertTrue(_testCase, colMap->find(dadSeq) == colMap->end() ||
                   colMap->find(dadSeq)->second->size() == 0);      
    }
    else
    {
      CuAssertTrue(_testCase, colMap->find(adamSeq) != colMap->end());
      ColumnIterator::DNASet* entry = colMap->find(adamSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      DNAIteratorConstPtr dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == adamSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i);

      CuAssertTrue(_testCase, colMap->find(grandpaSeq) != colMap->end());
      entry = colMap->find(grandpaSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == grandpaSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i - 4);

      CuAssertTrue(_testCase, colMap->find(dadSeq) != colMap->end());      
      entry = colMap->find(dadSeq)->second;
      CuAssertTrue(_testCase, entry->size() == 1);
      dna = entry->at(0);
      CuAssertTrue(_testCase, dna->getSequence() == dadSeq);
      CuAssertTrue(_testCase, dna->getArrayIndex() == i - 8);
    }

    colIterator->toRight();
  }
}

void ColumnIteratorPositionCacheTest::createCallBack(AlignmentPtr alignment)
{
}

void ColumnIteratorPositionCacheTest::checkCallBack(AlignmentConstPtr alignment)
{
  size_t trials = 10;
  size_t sizes[] = {10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 10000, 1000000};
  size_t entries = 10000;
  set<hal_index_t> truth;
  PositionCache cache;
  srand(time(NULL));

  for (size_t i = 0; i < trials; ++i)
  {
    for (size_t j = 0; j < entries; ++j)
    {
      hal_index_t val = (hal_index_t)rand() % sizes[i];
      bool r = truth.insert(val).second;
      bool r2 = cache.insert(val);
      CuAssertTrue(_testCase, r == r2);
      CuAssertTrue(_testCase, truth.size() == cache.size());
    }
    CuAssertTrue(_testCase, cache.check());
    for (size_t j = 0; j < entries * 2; ++j)
    {
      hal_index_t val = (hal_index_t)rand() % sizes[i];
      bool r = truth.find(val) != truth.end();
      bool r2 = cache.find(val);
      CuAssertTrue(_testCase, r == r2);
    }
    truth.clear();
    cache.clear();
  }
}

void halColumnIteratorBaseTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorBaseTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorDepthTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorDepthTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorDupTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorDupTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorInvTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorInvTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorGapTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorGapTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorMultiGapTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorMultiGapTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorMultiGapInvTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorMultiGapInvTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorPositionCacheTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorPositionCacheTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halColumnIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halColumnIteratorBaseTest);
  SUITE_ADD_TEST(suite, halColumnIteratorDepthTest); 
  SUITE_ADD_TEST(suite, halColumnIteratorDupTest);
  SUITE_ADD_TEST(suite, halColumnIteratorInvTest); 
  SUITE_ADD_TEST(suite, halColumnIteratorGapTest);
  SUITE_ADD_TEST(suite, halColumnIteratorMultiGapTest);
  SUITE_ADD_TEST(suite, halColumnIteratorMultiGapInvTest);
  SUITE_ADD_TEST(suite, halColumnIteratorPositionCacheTest);
  return suite;
}

