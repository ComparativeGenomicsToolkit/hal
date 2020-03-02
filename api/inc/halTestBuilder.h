/*
 * Copyright (C) 2012-2020 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * This module contains a set of class to simplify the building of test
 * alignments.  They let one specify the minimum about of information to
 * construct a test case in a textual manner and generate other data, such as
 * segments and unaligned genome.  This is designed to be easy for small
 * cases, not efficient.  Thus input consists of objects that are copied,
 * thus can be allocated on the stack.
 * 
 * It would be possible to have a text format file, such as JSON, as an
 * alternative to building in code.
 *
 * The *Spec objects are supplied by the client, they are extended as *Bld objects
 * with all the internal bookkeeping needed to build the HAL.
 */

#ifndef HALTESTBUILDER_H
#define HALTESTBUILDER_H
#include "halDefs.h"
#include <vector>
#include <map>

namespace hal {

    /* Class used to specify a genome */
    class GenomeSpec {
        public:
        GenomeSpec(const std::string& name,
                   const std::string& parentName="",
                   double branchLength=0.0):
            _name(name), _parentName(parentName), _branchLength(branchLength) {
        }
        const std::string& getName() const {
            return _name;
        }
        const std::string& getParentName() const {
            return _parentName;
        }
        double getBranchLength() const {
            return _branchLength;
        }
       
        private:
        const std::string _name;
        const std::string _parentName="";
        const double _branchLength;
    };
    typedef std::vector<GenomeSpec> GenomeSpecVec;

    /* Class to specify a sequence in a genome.  It does not specify the
     * bases, they are take from the alignment blocks or generated.
     * Note that these objects are always used by value, so they don't have
     * to be unique..*/
    class SeqSpec {
        public:
        SeqSpec(const std::string& genomeName,
                const std::string& name,
                hal_index_t length) :
            _genomeName(genomeName), _name(name), _fullName(genomeName + "." + name), _length(length) {
        }
        SeqSpec(const GenomeSpec& genome,
                const std::string& name,
                hal_index_t length) :
            _genomeName(genome.getName()), _name(name), _fullName(genome.getName() + "." + name), _length(length) {
        }
        const std::string& getGenomeName() const {
            return _genomeName;
        }
        const std::string& getName() const {
            return _name;
        }
        const std::string& getFullName() const {
            return _fullName;
        }
        hal_index_t getLength() const {
            return _length;
        }
        private:
        const std::string _genomeName;
        const std::string _name;
        const std::string _fullName;
        const hal_index_t _length;
    };
    typedef std::vector<SeqSpec> SeqSpecVec;
    
    /* Class to specify one component in a row block. */
    class RowSpec {
        public:
        /* Constructor.  The start and end coordinates are strand-specific,
         * as with MAF and must match the length of the bases in text. Specifying
         * both start and end is intended to be more obvious, even though end is redundant.
         * The text string contains bases, N, or '-' characters, */
        RowSpec(const SeqSpec& seq, hal_index_t start, hal_index_t end, char strand, 
                const std::string& columns):
            _genomeName(seq.getGenomeName()), _name(seq.getName()), _fullName(seq.getFullName()),
            _start(start), _end(end), _strand(strand), _columns(columns) {
        }
        std::string getDesc() const {
            return _fullName + ":" + std::to_string(_start) + "-" + std::to_string(_end) + "/" + std::string(1, _strand);
        }
        const std::string& getGenomeName() const {
            return _genomeName;
        }
        const std::string& getName() const {
            return _name;
        }
        const std::string& getFullName() const {
            return _fullName;
        }
        hal_index_t getStart() const {
            return _start;
        }
        hal_index_t getEnd() const {
            return _end;
        }
        char getStrand() const {
            return _strand;
        }
        const std::string getColumns() const {
            return _columns;
        }

        private:
        const std::string _genomeName;
        const std::string _name;
        const std::string _fullName;
        const hal_index_t _start;
        const hal_index_t _end;
        const char _strand;
        const std::string _columns;
    };
    typedef std::vector<RowSpec> RowSpecVec;


    /* A block of the alignment, which is a vector of RowSpec objects */
    class BlockSpec: public RowSpecVec {
        public:
        BlockSpec(const std::initializer_list<RowSpec>& src) {
            for (auto it = src.begin(); it != src.end(); it++) {
                push_back(*it);
            }
        }

    };
    typedef std::vector<BlockSpec> BlockSpecVec;

    /** 
     * Class used to collect the specifications and then validate and generate
     * the alignment.
     **/
    class TestBuilder {
        public:

        /* parent must be define before children */
        void addGenome(const GenomeSpec& genome);

        /* parent must be define before children */
        void addGenomes(const GenomeSpecVec& genomes);

        /* genome must already be defined */
        void addSeq(const SeqSpec& seq);

        /* genome must already be defined */
        void addSeqs(const SeqSpecVec& seqs);
        
        /* sequence must already be defined */
        void addBlock(const BlockSpec& block);

        /* sequence must already be defined */
        void addBlocks(const BlockSpecVec& blocks);

        void add(const GenomeSpecVec& genomes,
                 const SeqSpecVec& seqs,
                 const BlockSpecVec& blocks = BlockSpecVec()) {
            addGenomes(genomes);
            addSeqs(seqs);
            addBlocks(blocks);
        }
        
        /* build the alignment */
        void build(Alignment *alignment);
        
        /* print for debugging */
        void print(std::ostream& os) const;
        
        private:
        // forwards
        class GenomeBld;
        class SeqBld;
        class RowBld;
        class BlockBld;
        class SeqSegBld;
        class SegBld;

        typedef std::vector <GenomeBld*> GenomeBldVec;
        typedef std::vector<SeqBld*> SeqBldVec;
        
        GenomeBld* findGenomeBld(const std::string& name) {
            auto it = _genomeBlds.find(name);
            return (it != _genomeBlds.end()) ? it->second : nullptr;
        }
        GenomeBld* getGenomeBld(const std::string& name);

        SeqBld* findSeqBld(const std::string& fullName) {
            auto it = _seqBlds.find(fullName);
            return (it != _seqBlds.end()) ? it->second : nullptr;
        }

        SeqBld* getSeqBld(const std::string& fullName);

        SeqBld* checkRow(const RowSpec& row);
        void fillBases(RowBld* rowBld);
        void fillBasesFromBlocks();
        void fakeUnfilledBases();
        void buildSeqs();
        void buildSegments();

        /* pointers are stored so don't worry about objects being relocated */
        std::map<const std::string, GenomeBld*> _genomeBlds;
        std::map<const std::string, SeqBld*> _seqBlds;
        std::vector<BlockBld*> _blockBlds;
    };
    
}
#endif
