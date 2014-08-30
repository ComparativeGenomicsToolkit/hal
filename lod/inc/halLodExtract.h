/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODEXTRACT_H
#define _HALLODEXTRACT_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"
#include "halLodGraph.h"

namespace hal {

/** This is the class that is responsible for reading and writing the 
 * HAL alignments: for each internal node to be written, it builds 
 * a sequence graph from the input alignment, then writes that sequence
 * graph as an internal node (and immediate descendants) in the output
 * alignment.  We always make the assumption that we are starting from
 * scratch (ie the output alignment does not already exist). 
 *
 * The output alignment is created from an arbitrary subset of genomes from
 * the input, linked together in an arbitrary tree.  By default, the 
 * identical tree is used. */
class LodExtract
{
public:
   
   LodExtract();
   ~LodExtract();

   void createInterpolatedAlignment(AlignmentConstPtr inAlignment,
                                    AlignmentPtr outAlignment,
                                    double scale,
                                    const std::string& tree,
                                    const std::string& rootName,
                                    bool keepSequences,
                                    bool allSequences,
                                    double probeFrac,
                                    double minSeqFrac);
   
   
protected:

   typedef std::set<const LodSegment*, LodSegmentPLess> SegmentSet;
   typedef std::map<const Genome*, SegmentSet*> SegmentMap;

protected:

   void createTree(const std::string& tree, const std::string& rootName);
   void convertInternalNode(const std::string& genomeName, double scale);
   void countSegmentsInGraph(
     std::map<const Sequence*, hal_size_t>& segmentCounts);
   void writeDimensions(
     const std::map<const Sequence*, hal_size_t>& segmentCounts, 
     const std::string& parentName,
     const std::vector<std::string>& childNames);
   void writeSequences(const Genome* inParent, 
                      const std::vector<const Genome*>& inChildren);
   void writeSegments(const Genome* inParent, 
                      const std::vector<const Genome*>& inChildren);
   void writeUnsampledSequence(const Sequence* outSequence,
                               SegmentIteratorPtr outSegment);
   void writeHomologies(const Genome* inParent, 
                        const std::vector<const Genome*>& inChildren);
   void updateBlockEdges(const Genome* inParentGenome,
                         SegmentMap& segMap,
                         const LodBlock* block,
                         BottomSegmentIteratorPtr bottom,
                         TopSegmentIteratorPtr top);
   void writeParseInfo(Genome* genome);
   hal_size_t getMinAvgBlockSize(
     const Genome* inParent,
     const std::vector<const Genome*>& inChildren,
     const Genome* inGrandParent) const;

   
   AlignmentConstPtr _inAlignment;
   AlignmentPtr _outAlignment;

   LodGraph _graph;
   bool _keepSequences;
   bool _allSequences;
   double _probeFrac;
   double _minSeqFrac;
};

}

#endif
