/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include "hal2psl.h"
#include "psl.h"


void hal::Hal2Psl::storePslResults(std::vector<PslBlock>& pslBlocks){
    
    BedList::iterator i = _outBedLines.begin();
    for (; i != _outBedLines.end(); ++i)
    {
      if (_addExtraColumns == false)
      {
        i->_extra.clear();
      }
      makeUpPsl(i->_psl, i->_blocks, i->_strand, i->_start, i->_chrName, 
                 pslBlocks);
    }
    
}
std::vector<PslBlock> hal::Hal2Psl::convert2psl(hal::AlignmentConstPtr alignment,
                       const hal::Genome* srcGenome,
                       const hal::Genome* tgtGenome,
                       const std::string srcChrom){
    
    std::vector<PslBlock> pslBlocks;
    if (srcGenome->getNumSequences() > 0){
        _srcGenome = srcGenome;
        _tgtGenome = tgtGenome; 
        _coalescenceLimit = NULL;
        _traverseDupes = true;
        _addExtraColumns = false;
        _missedSet.clear();
        _tgtSet.clear();
        _tgtSet.insert(tgtGenome);
        hal::SequenceIteratorConstPtr seqIt = srcGenome->getSequenceIterator();
        hal::SequenceIteratorConstPtr seqEnd = srcGenome->getSequenceEndIterator();
        for (; !seqIt->equals(seqEnd); seqIt->toNext()) {
          //std::cout << "aa" << srcChrom << std::endl;
          _outBedLines.clear();
          _srcSequence = seqIt->getSequence(); 
          auto chrName = _srcSequence->getName(); 
          if (srcChrom != "\"\""  && srcChrom != chrName) {continue;}
          //else {std::cout << srcChrom << std::endl;}
          _bedLine._start = 0;
          _bedLine._end = _srcSequence->getSequenceLength(); 
          _mappedBlocks.clear();
          _outPSL = true;
          visitBegin();
          liftInterval(_mappedBlocks);  
          if (_mappedBlocks.size()) {
            assignBlocksToIntervals();
          }
          //_outBedLines.sort(BedLineSrcLess());
          storePslResults(pslBlocks);
          cleanResults();
        }
    }
    return pslBlocks;
}


void hal::Hal2Psl::makeUpPsl(const std::vector<hal::PSLInfo>& vpsl,
                                const std::vector<hal::BedBlock>& blocks, 
                                const char strand,
                                const hal_index_t start,
                                const std::string chrName,
                                std::vector<PslBlock>& pslBlocks){
  assert(vpsl.size() == 1);
  const hal::PSLInfo& psl = vpsl[0];

  for (size_t i = 0; i < blocks.size(); ++i) {
    PslBlock b = PslBlock();
    b.size = blocks[i]._length;
    b.qName = psl._qSeqName;
    b.qSize = psl._qSeqSize;
    
    b.qStart = psl._qBlockStarts[i] - psl._qChromOffset;
    b.qEnd = b.qStart + b.size;
    if (psl._qStrand == '-'){
        auto posStart = b.qStart;
        b.qStart = psl._qSeqSize - posStart - b.size;
        b.qEnd = psl._qSeqSize - posStart;
    }
    
    b.tStart = blocks[i]._start + start;
    b.tEnd = b.tStart + b.size;
    if (strand == '-')
    {
      auto posStart = b.tStart;
      b.tStart = psl._tSeqSize - posStart - b.size;
      b.tEnd = psl._tSeqSize - posStart;
    }
    std::stringstream ss;
    ss << psl._qStrand << strand;
    ss >> b.strand;
    //b.strand = psl._qStrand+'/'+strand;
    b.tSize = psl._tSeqSize;
    b.qSize = psl._qSeqSize;
    b.tName = chrName;
    pslBlocks.push_back(b);
  }
}

