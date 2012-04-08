/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENT_H
#define _HALALIGNMENT_H

#include <string>
#include "halDefs.h"

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

   /** Close the alignment */
   virtual void close() = 0;
   
   /** Add a new genome to the alignment
    * @param name of new genome in alignment (must be unique)
    * @param path subdirectory of new genome (empty for default)
    * @param parentName name of parent genome (empty if adding root)
    * @param childNames paths of child genomes (if any) */
   virtual GenomePtr addGenome(const std::string& name,
                               const std::string& path,
                               const std::string& parentName,
                               const std::vector<std::string>& childNames) = 0;

   /** Remove a genome from the alignment 
    * @param path Path of genome to remove */
   virtual void removeGenome(const std::string& name) = 0;

   /** Open an existing genome for reading and updating
    * @path Path of genome to open */
   virtual GenomeConstPtr openConstGenome(const std::string& path) const = 0;

   /** Open an exsting genome for reading
    * @path Path of genome to open */
   virtual GenomePtr openGenome(const std::string& path) = 0;

   /** Get name of root genome (empty string for empty alignment) */
   virtual std::string getRootName() const = 0;

   /** Get name of parent genome in the phylogeny (empty string for root)
    * @param name Name of genome */
   virtual std::string getParentName(const std::string& name) const = 0;
   
   /** Get names of child genomes in the phylogeny (empty vector for leaves)
    * @param name Name of genome */
   virtual std::vector<std::string> 
   getChildNames(const std::string& name) const= 0;

   /** Get the names of all leaves below a given genome 
    * @param name Name of genome */
   virtual std::vector<std::string>
   getLeafNamesBelow(const std::string& name) const = 0;

   /** Get Alignment's metadata */
   virtual MetaDataPtr getMetaData() = 0;

   /** Get read-only instance of Alignment's metadata */
   virtual MetaDataConstPtr getMetaData() const = 0;

   /** Get path of a genome 
    * @param name Name of genome whose path we want */
   virtual const std::string& getPath(const std::string& name) const = 0;
};

inline Alignment::~Alignment() {}
}
#endif
