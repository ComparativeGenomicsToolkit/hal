/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENT_H
#define _HALALIGNMENT_H

#include "halDefs.h"
#include <string>
#include <vector>

namespace hal {
    /* thrown when genome not found in alignment */
    class GenomeNotFoundException: public hal_exception {
        public:
        GenomeNotFoundException(const std::string &name):
            hal_exception("Genome not found: '" + name + "'") {
        }
    };

    /**
     * Interface for a hierarhcical alignment.  Responsible for creating
     * and accessing genomes and tree information.  Accesssing a HAL file must
     * always start with an Alignment object
     */
    class Alignment {
      public:
        /** Destructor */
        virtual ~Alignment() {
        }

        /** Close the alignment */
        virtual void close() = 0;

        /** Close the alignment */
        void close() const {
            const_cast<Alignment*>(this)->close();
        }

        /** Get the name of the storage engine.
         **/
        virtual const std::string &getStorageFormat() const = 0;

        /** Add a new genome to the alignment
         * @param name name of new genome in alignment (must be unique)
         * @param parent name of parent genome in tree (must exist)
         * @param branchLength distance between new genome and parent */
        virtual Genome *addLeafGenome(const std::string &name, const std::string &parentName, double branchLength) = 0;

        /** Add a new genome as root to the alignment.  The previous root
         * (if exists) will be a child of the new genome
         * @param name name of new genome in alignment (must be unique)
         * @param branchLength distance between new genome and previous root
         * (if exists)*/
        virtual Genome *addRootGenome(const std::string &name, double branchLength = 0) = 0;

        /** Remove a leaf genome from the alignment
         * @param path Path of genome to remove */
        virtual void removeGenome(const std::string &name) = 0;

        /** Insert a new genome between a node and its child
         * @param name Name of new genome
         * @param parentName Name of the existing parent genome
         * @param child Name of the existing child genome
         * @param upperBranchLength Length of parent-insert branch
         * (length of insert-child branch will be inferred) */
        virtual Genome *insertGenome(const std::string &name, const std::string &parentName, const std::string &childName,
                                     double upperBranchLength) = 0;

        /** Open an exsting genome for reading
         * @param name Name of genome to open
         * @return Genome object or NULL is not found
         * DO NOT DELETE THE RETURNED OBJECT, USE closeGenome() */
        virtual const Genome *openGenome(const std::string &name) const = 0;

        /** Open an existing genome for reading and updating
         * @param name Name of genome to open.
         * @return Genome object or NULL is not found
         * DO NOT DELETE THE RETURNED OBJECT, USE closeGenome() */
        virtual Genome *openGenome(const std::string &name) = 0;

        /** Open an exsting genome for reading, checking for missing genome
         * @param name Name of genome to open
         * @return Genome object 
         * @throws GenomeNotFoundException
         * DO NOT DELETE THE RETURNED OBJECT, USE closeGenome() */
        const Genome *openGenomeCheck(const std::string &name) const {
            const Genome *genome = openGenome(name);
            if (genome == NULL) {
                throw GenomeNotFoundException(name);
            }
            return genome;
        }

        /** Open an existing genome for reading and updating, checking for missing genome
         * @param name Name of genome to open.
         * @return Genome object 
         * @throws GenomeNotFoundException
         * DO NOT DELETE THE RETURNED OBJECT, USE closeGenome() */
        Genome *openGenomeCheck(const std::string &name) {
            Genome *genome = openGenome(name);
            if (genome == NULL) {
                throw GenomeNotFoundException(name);
            }
            return genome;
        }

        /** Close an open genome.  All pointers to this genome become
         * invalid and openGenome needs to be called again to access it
         * @param genome Genome to close */
        virtual void closeGenome(const Genome *genome) const = 0;

        /** Get name of root genome (empty string for empty alignment) */
        virtual std::string getRootName() const = 0;

        /** Get name of parent genome in the phylogeny (empty string for root)
         * @param name Name of genome */
        virtual std::string getParentName(const std::string &name) const = 0;

        /** Get the branch length between two genomes in the phylogeny
         * @param parentName name of parent genome
         * @param childName name of child genome */
        virtual double getBranchLength(const std::string &parentName, const std::string &childName) const = 0;

        /** Change the branch length between two genomes in the phylogeny
         * @param parentName name of parent genome
         * @param childName name of child genome */
        virtual void updateBranchLength(const std::string &parentName, const std::string &childName, double length) = 0;

        /** Get names of child genomes in the phylogeny (empty vector for leaves)
         * @param name Name of genome */
        virtual std::vector<std::string> getChildNames(const std::string &name) const = 0;

        /** Get the names of all leaves below a given genome
         * @param name Name of genome */
        virtual std::vector<std::string> getLeafNamesBelow(const std::string &name) const = 0;

        /** Get the number of genomes (including internal ancestors)
         * in the alignment */
        virtual hal_size_t getNumGenomes() const = 0;

        /** Get Alignment's metadata */
        virtual MetaData *getMetaData() = 0;

        /** Get read-only instance of Alignment's metadata */
        virtual const MetaData *getMetaData() const = 0;

        /** Get a newick-formatted phylogeny to the alignment */
        virtual std::string getNewickTree() const = 0;

        /** Get version used to create the file */
        virtual std::string getVersion() const = 0;

        /** Is this file open for read-only? */
        virtual bool isReadOnly() const = 0;

        /** Replace the newick tree with a new string */
        virtual void replaceNewickTree(const std::string &newick) = 0;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
