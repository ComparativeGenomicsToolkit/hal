/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENT_H
#define _HALALIGNMENT_H

#include <string>
#include "halDefs.h"
#include "halMetaData.h"
#include "halGenome.h"

namespace hal {

/** 
 * Interface for a hierarhcical alignment
 * Right now genomes are assigned paths within the alignment
 * in order to allow for an arbitrary grouping structure.  This
 * may be refined at some point. 
 * TODO: revise tree-related stuff (maybe move to genome class)
 */
class Alignment
{
public:

   /** Destructor */
   virtual ~Alignment() = 0;

   /** Create new alignment (overwriting any existing alignments)
    * @param alignmentPath location on disk */
   virtual void createNew(const std::string& alignmentPath) = 0;

   /** Open an existing alignment 
    * @param alignmentPath location on disk 
    * @param readOnly access flag */
   virtual void open(const std::string& alignmentPath, 
                     bool readOnly) = 0;
   
   /** Add a new genome to the alignment
    * @param path location of new genome in alignment
    * @param parentPath path of parent genome (null if adding root)
    * @param childPaths paths of child genomes (if any) */
   virtual GenomePtr addGenome(const std::string& path, 
                               const std::string* parentPath,
                               const std::vector<std::string>& childPaths) = 0;

   /** Remove a genome from the alignment 
    * @param path Path of genome to remove */
   virtual void removeGenome(const std::string& path) = 0;

   /** Open an existing genome for reading and updating
    * @path Path of genome to open */
   virtual GenomeConstPtr openConstGenome(const std::string& path) const = 0;

   /** Open an exsting genome for reading
    * @path Path of genome to open */
   virtual GenomePtr openGenome(const std::string& path) = 0;

   /** Get path of parent genome in the phylogeny (empty string for root)
    * @param path Path of genome */
   virtual std::string getParent(const std::string& path) = 0;
   
   /** Get paths of child genomes in the phylogeny (empty vector for leaves)
    * @param path Path of genome */
   virtual std::vector<std::string> getChildren(const std::string& path) = 0;

   /** Get Alignment's metadata */
   virtual MetaDataPtr getMetaData() = 0;

   /** Get read-only instance of Alignment's metadata */
   virtual MetaDataConstPtr getMetaData() const = 0;
};


}
#endif
