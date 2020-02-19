/*
 * Copyright (C) 2012-2020 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * Module to simplify the building of test alignments.
 */
#include "halTestBuilder.h"
#include "halCommon.h"
#include <ostream>
using namespace std;
using namespace hal;

/* GenomeSpec with state */
class TestBuilder::GenomeBld: public GenomeSpec {
    public:
    GenomeBld(const GenomeSpec& genome):
        GenomeSpec(genome) {
    }

    void addSeqBld(SeqBld* seqBld) {
        _seqBlds.push_back(seqBld);
    }
    void print(ostream& os) const;
    
    private:
    vector<SeqBld*> _seqBlds;
};

/* SeqSpec with state */
class TestBuilder::SeqBld: public SeqSpec {
    public:
    SeqBld(const SeqSpec& seq,
           GenomeBld* genomeBld):
        SeqSpec(seq), _genomeBld(genomeBld) {
        _bases.resize(seq.getLength());
    }
    string* getBases() {
        return &_bases;
    }
    
    void print(ostream& os) const {
        os << "seq: " << getFullName() << " len: " << getLength() << " " << _bases << endl;
    }
    
    private:
    GenomeBld* _genomeBld;
    string _bases;
};

/* RowSpec with state */
class TestBuilder::RowBld: public RowSpec {
    public:
    RowBld(const RowSpec& row,
           SeqBld* seqBld):
        RowSpec(row), _seqBld(seqBld) {
    }
    SeqBld* getSeqBld() {
        return _seqBld;
    }
    
    private:
    SeqBld* _seqBld;
};

/* BlockSpec with state */
class TestBuilder::BlockBld: public vector<RowBld*> {
    public:
    BlockBld(vector<RowBld*> rowBlds):
        vector<RowBld*>(rowBlds) {
    }
    private:
};


void TestBuilder::GenomeBld::print(ostream& os) const {
    os << "gnome: " << getName() << "  parent: " << getParentName() << endl;
    for (auto it = _seqBlds.begin(); it != _seqBlds.end(); it++) {
        os << "    ";
        (*it)->print(os);
        }
}


static hal_index_t countBases(const string& s) {
    hal_index_t cnt = 0;
    for (auto i = 0; i < s.size(); i++) {
        if (s[i] != '-') {
            cnt++;
        }
    }
    return cnt;
}

void TestBuilder::addGenome(const GenomeSpec& genome) {
    if (_genomeBlds.find(genome.getName()) != _genomeBlds.end()) {
        throw hal_exception("test genome already defined: " + genome.getName());
    }
    _genomeBlds.insert({genome.getName(), new GenomeBld(genome)});
}
    
void TestBuilder::addGenomes(const GenomeSpecVec& genomes) {
    for (auto it = genomes.begin(); it != genomes.end(); it++) {
        addGenome(*it);
    }
}

void TestBuilder::addSeq(const SeqSpec& seq) {
    if (_seqBlds.find(seq.getFullName()) != _seqBlds.end()) {
        throw hal_exception("sequence already defined: " + seq.getFullName());
    }
    GenomeBld* genomeBld = findGenomeBld(seq.getGenomeName());
    if (genomeBld == nullptr) {
        throw hal_exception("genome not defined for sequence: " + seq.getFullName());
    }
    SeqBld* seqBld = new SeqBld(seq, genomeBld);
    genomeBld->addSeqBld(seqBld);
    _seqBlds.insert({seq.getFullName(), seqBld});
    
}

void TestBuilder::addSeqs(const SeqSpecVec& seqs) {
    for (auto it = seqs.begin(); it != seqs.end(); it++) {
        addSeq(*it);
    }
}

TestBuilder::SeqBld* TestBuilder::checkRow(const RowSpec& row) {
    SeqBld* seqBld = findSeqBld(row.getFullName());
    if (seqBld == nullptr) {
        throw hal_exception("sequence in block row not define: " + row.getFullName());
    }
    if (row.getStart() >= row.getEnd()) {
        throw hal_exception("negative or zero length block row: " + row.getFullName());
    }
    if (row.getEnd() > seqBld->getLength()) {
        throw hal_exception("block row out of sequence bound: " + row.getFullName());
    }
    char strand = row.getStrand();
    if (not ((strand == '+') or (strand == '-'))) {
        throw hal_exception("block row has invalid strand: " + row.getFullName());
    }

    hal_index_t baseCnt = countBases(row.getColumns());
    if (baseCnt != (row.getEnd() - row.getStart())) {
        throw hal_exception("block row has length " + to_string(row.getEnd() - row.getStart())
                            + " however it has a base count of " + to_string(baseCnt));
    }
    return seqBld;
}

void TestBuilder::addBlock(const BlockSpec& block) {
    vector<RowBld*> rowBlds;
    for (auto it = block.begin(); it != block.end(); it++) {
        SeqBld* seqBld = checkRow(*it);
        rowBlds.push_back(new RowBld(*it, seqBld));
    }
    _blockBlds.push_back(new BlockBld(rowBlds));
}

void TestBuilder::addBlocks(const BlockSpecVec& blocks) {
    for (auto it = blocks.begin(); it != blocks.end(); it++) {
        addBlock(*it);
    }
}

void TestBuilder::fillBases(RowBld* rowBld) {
    SeqBld *seqBld = rowBld->getSeqBld();
    string columns(rowBld->getColumns());
    hal_index_t iSeq, dir;
    if (rowBld->getStrand() == '+') {
        iSeq = rowBld->getStart();
        dir = 1;
    } else {
        complement(columns);
        iSeq = (seqBld->getLength() - rowBld->getStart()) - 1;
        dir = -1;
    }
    string& bases = *seqBld->getBases();
    for (hal_index_t iCol = 0; iCol < columns.size(); iCol++) {
        if (columns[iCol] != '-') {
            if ((bases[iSeq] != '\0') and (bases[iSeq] != columns[iCol])) {
                throw hal_exception("block sequence " + rowBld->getDesc() + " sets base to '" + string(1, bases[iSeq])
                                    + "', was set by another block to '" + string(1, columns[iCol]) + "'");
            }
            bases[iSeq] = columns[iCol];
            iSeq += dir;
        }
    }
}

void TestBuilder::fillBasesFromBlocks() {
    for (auto bit = _blockBlds.begin(); bit != _blockBlds.end(); bit++) {
        for (auto rit = (*bit)->begin(); rit != (*bit)->end(); rit++) {
            fillBases(*rit);
        }
    }
}

// set bases not in any of the blocks to N
void TestBuilder::fakeUnfilledBases() {
    for (auto it = _seqBlds.begin(); it != _seqBlds.end(); it++) {
        string& bases = *(it->second->getBases());
        for (hal_index_t i = 0; i < bases.size(); i++) {
            if (bases[i] == '\0') {
                bases[i] = 'N';
            }
        }
    }
}

void TestBuilder::buildSeqs() {
    fillBasesFromBlocks();
    fakeUnfilledBases();
}

void TestBuilder::build(Alignment *alignment) {
    buildSeqs();
}

void TestBuilder::print(ostream& os) const {
    for (auto it = _genomeBlds.begin(); it != _genomeBlds.end(); it++) {
        it->second->print(os);
    }
}
