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
                                    hal_size_t step,
                                    const std::string& tree);
   
   
protected:

   void createTree(const std::string& tree);
   void convertInternalNode(const std::string& genomeName, hal_size_t step);
   void countSegmentsInGraph(
     std::map<const Sequence*, hal_size_t>& segmentCounts);
   void writeDimensions(
     const std::map<const Sequence*, hal_size_t>& segmentCounts, 
     const std::string& parentName,
     const std::vector<std::string>& childNames);

   
   AlignmentConstPtr _inAlignment;
   AlignmentPtr _outAlignment;

   LodGraph _graph;
};

}

#endif
