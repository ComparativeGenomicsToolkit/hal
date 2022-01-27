/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hal.h"
#include "halCLParser.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unordered_set>

using namespace std;
using namespace hal;

static void genome2PAF(ostream& outStream, const Genome* genome, bool fullNames);

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("inHalPath", "input hal file");
    optionsParser.addOption("rootGenome", 
                            "process only genomes in clade with specified root"
                            " (HAL root if empty)", 
                            "\"\"");
    optionsParser.addOptionFlag("onlySequenceNames", "use only sequence names "
                                "for output names.  By default, the UCSC convention of Genome.Sequence "
                                "is used",
                                false);
    optionsParser.setDescription("Export pairwise alignment (with no softclips) of each branch to PAF");
}

int main(int argc, char **argv) {
    CLParser optionsParser;
    initParser(optionsParser);
    
    string halPath;
    string rootGenomeName;
    bool fullNames;

    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("inHalPath");
        rootGenomeName = optionsParser.getOption<string>("rootGenome");
        fullNames = !optionsParser.getFlag("onlySequenceNames");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    try {
        AlignmentConstPtr alignment(openHalAlignment(halPath, &optionsParser));
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("input hal alignmenet is empty");
        }
        
        // root is specified either by the parameter or as the alignment root
        // by default
        const Genome* rootGenome = NULL;
        if (rootGenomeName != "\"\"") {
            rootGenome = alignment->openGenome(rootGenomeName);
        } else {
            rootGenome = alignment->openGenome(alignment->getRootName());
        }
        if (rootGenome == NULL) {
            throw hal_exception(string("Root genome, ") + rootGenomeName + 
                                ", not found in alignment");
        }
        const Genome* parentGenome = rootGenome;

        vector<string> childs = alignment->getChildNames(rootGenome->getName());
        deque<string> queue(childs.begin(), childs.end());

        while (!queue.empty()) {
            string childName = queue.front();
            queue.pop_front();
            const Genome* childGenome = alignment->openGenome(childName);
            string parentName = alignment->getParentName(childName);
            if (parentName != parentGenome->getName()) {
                alignment->closeGenome(parentGenome);
                parentGenome = childGenome->getParent();
            }

            genome2PAF(cout, childGenome, fullNames);

            childs = alignment->getChildNames(childName);
            for (int i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }

            alignment->closeGenome(childGenome);
        }
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
     
    return 0;
}

/// scan to next match, returning false if not found
static bool nextMatch(const TopSegmentIteratorPtr& topIt1, const BottomSegmentIteratorPtr& botIt1,
                      TopSegmentIteratorPtr& topIt2,  BottomSegmentIteratorPtr& botIt2) {
    // set the second iterators to match the first, without re-allocating
    topIt2->copy(topIt1);
    botIt2->copy(botIt2);

    // scan til next match
    for (topIt2->toRight(); not topIt2->atEnd(); topIt2->toRight()) {
        if (topIt2->tseg()->hasParent()) {
            // return on any kind of match
            botIt2->toParent(topIt2);
            return true;
        }
    }

    topIt2->copy(topIt1);
    botIt2->copy(botIt1);

    return false;
}

/// check if given intervals form a contiguous block (m), insertion (i), deletion (d) or anything else (o)
static char blockCat(const TopSegmentIteratorPtr& topIt1, const BottomSegmentIteratorPtr& botIt1,
                     const TopSegmentIteratorPtr& topIt2, const BottomSegmentIteratorPtr& botIt2,
                     TopSegmentIteratorPtr& topIt3, BottomSegmentIteratorPtr& botIt3,
                     const unordered_set<hal_index_t>& parentSet) {
    // i: insertion d: deletion m: contiguous match o: other

    if (topIt1->getSequence() == topIt2->getSequence() &&
        botIt1->getSequence() == botIt2->getSequence() &&
        botIt1->getReversed() == botIt2->getReversed()) {

        assert(topIt2->getArrayIndex() > topIt1->getArrayIndex());

        bool topAdj = topIt2->getArrayIndex() == topIt1->getArrayIndex() + 1;
        bool botAdj = (botIt1->getReversed() && botIt1->getArrayIndex() == botIt2->getArrayIndex() + 1) ||
            (!botIt1->getReversed() && botIt2->getArrayIndex() == botIt1->getArrayIndex() + 1);
        
        // check match
        if (topAdj && botAdj) {
            return 'm';
        }

        // check insertion
        if (topIt1->getArrayIndex() + 1 < topIt2->getArrayIndex() && botAdj) {
            bool is_insert = true;
            for (topIt3->copy(topIt1); !topIt3->equals(topIt2) && is_insert; topIt3->toRight()) {
                if (!topIt3->equals(topIt1) && !topIt3->equals(topIt2)) {
                    is_insert = !topIt3->tseg()->hasParent();
                }
            }
            if (is_insert) {
                return 'i';
            }
        }

        // check deletion
        if (topAdj &&
            ((botIt1->getReversed() && botIt1->getArrayIndex() > botIt2->getArrayIndex() + 1) ||
             (!botIt1->getReversed() && botIt1->getArrayIndex() + 1 < botIt2->getArrayIndex()))) {
            bool is_del = true;
            for (botIt3->copy(botIt1); !botIt3->equals(botIt2) && is_del; botIt3->toRight()) {
                if (!botIt3->equals(botIt1) && !botIt3->equals(botIt2)) {
                    is_del = !botIt3->bseg()->hasChildG(topIt1->getGenome()) && !parentSet.count(botIt3->bseg()->getArrayIndex());
                }
            }
            if (is_del) {
                return 'd';
            }
        } 
    }
    
    return 'o';
}

static size_t countSnps(TopSegmentIteratorPtr& topIt, BottomSegmentIteratorPtr& botIt) {
    string top_bases;
    topIt->getString(top_bases);
    string bot_bases;
    botIt->getString(bot_bases);
    assert(top_bases.length() == bot_bases.length());
    size_t num_snps = 0;
    for (size_t i = 0; i < top_bases.length(); ++i) {
        if (fastUpper(top_bases[i]) != fastUpper(bot_bases[i])) {
            ++num_snps;
        }
    }
    return num_snps;
}

void genome2PAF(ostream& outStream, const Genome* genome, bool fullNames) {
    TopSegmentIteratorPtr topIt1 = genome->getTopSegmentIterator();
    TopSegmentIteratorPtr topIt2 = genome->getTopSegmentIterator();
    TopSegmentIteratorPtr topIt3 = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIt1 = genome->getParent()->getBottomSegmentIterator();
    BottomSegmentIteratorPtr botIt2 = genome->getParent()->getBottomSegmentIterator();
    BottomSegmentIteratorPtr botIt3 = genome->getParent()->getBottomSegmentIterator();
    bool found_match = false;

    // this is used to remember which bottom segments have children in presence of duplications
    unordered_set<hal_index_t> parentSet;
    for (; not topIt1->atEnd(); topIt1->toRight()) {
        if (topIt1->tseg()->hasNextParalogy() && topIt1->tseg()->isCanonicalParalog()) {
            parentSet.insert(topIt1->tseg()->getParentIndex());
        }
    }
    topIt1->copy(topIt2);

    // find the first match
    if (topIt1->tseg()->hasParent()) {
        // it's either at position 0
        botIt1->toParent(topIt1);
        botIt2->toParent(topIt2);
        topIt2->copy(topIt1);
        botIt2->copy(botIt1);
        found_match = true;
    }
    if (!found_match) {
        // or we have to seek for it
        found_match = nextMatch(topIt1, botIt1, topIt2, botIt2);
        if (!found_match) {
            cerr << "Warning [hal2paf]: no alignment blocks found for genome " << genome->getName() << endl;
            // don't bother printing out empty records
            return;
        }
        topIt1->copy(topIt2);
        botIt1->copy(botIt2);
    }

    // start up a fresh paf line
    string queryName = fullNames ? topIt1->getSequence()->getFullName() : topIt1->getSequence()->getName();
    size_t queryLength = topIt1->getSequence()->getSequenceLength();
    hal_index_t queryStart = topIt1->getStartPosition();
    hal_index_t queryEnd = topIt1->getEndPosition();
    string targetName = fullNames ? botIt1->getSequence()->getFullName() : botIt1->getSequence()->getName();
    size_t targetLength = botIt1->getSequence()->getSequenceLength();
    hal_index_t targetStart = botIt1->bseg()->getStartPosition();
    hal_index_t targetEnd = botIt1->bseg()->getEndPosition();
    size_t matches = topIt1->getLength();
    size_t snps = countSnps(topIt1, botIt1);
    size_t gaps = 0;
    size_t runningMatch = topIt1->getLength();
    vector<string> cigar;
        
    // go forward
    while (!cigar.empty() || runningMatch > 0) {

#ifdef debug
        cerr << "t1 " << *topIt1 << endl;
        cerr << "b1 " << *botIt1 << endl;
        cerr << "t2 " << *topIt2 << endl;
        cerr << "b2 " << *botIt2 << endl;        
#endif
        found_match = nextMatch(topIt1, botIt1, topIt2, botIt2);
        bool newLine = !found_match;
        char cat = 'o';
        if (!newLine) {
            cat = blockCat(topIt1, botIt1, topIt2, botIt2, topIt3, botIt3, parentSet);
            newLine = cat == 'o';
        }
#ifdef debug
        cerr << "after match with cat " << string(1, cat) << " and found match " << found_match << endl;
        cerr << "t1 " << *topIt1 << endl;
        cerr << "b1 " << *botIt1 << endl;
        cerr << "t2 " << *topIt2 << endl;        
        cerr << "b2 " << *botIt2 << endl;
        cerr << "---" << endl;
#endif
            
        if (newLine) {
            // resolve the running match
            if (runningMatch > 0) {
                cigar.push_back(std::to_string(runningMatch) + "M");
                runningMatch = 0;
            }

            // make our cigar
            string cigar_string;
            if (botIt1->getReversed()) {
                for (vector<string>::reverse_iterator ci = cigar.rbegin(); ci != cigar.rend(); ++ci) {
                    cigar_string += *ci;
                }
            } else {
                for (vector<string>::iterator ci = cigar.begin(); ci != cigar.end(); ++ci) {
                    cigar_string += *ci;
                }
            }
            
            // write out the current paf line
            outStream << queryName << "\t" << queryLength << "\t"
                      << (queryStart - topIt1->getSequence()->getStartPosition()) << "\t"
                      << (queryEnd - topIt1->getSequence()->getStartPosition() + 1) << "\t"
                      << (botIt1->getReversed() ? "-" : "+") << "\t"
                      << targetName << "\t" << targetLength << "\t"
                      << (targetStart - botIt1->getSequence()->getStartPosition()) << "\t"
                      << (targetEnd - botIt1->getSequence()->getStartPosition() + 1) << "\t"
                      << (matches - snps)<< "\t"
                      << (matches + gaps) << "\t"
                      << 255 << "\t"
                      << "cg:Z:" << cigar_string << endl;

            cigar.clear();
            
            if (found_match) {
                // start a new cigar
                runningMatch = topIt2->getLength();
                queryName = fullNames ? topIt2->getSequence()->getFullName() : topIt2->getSequence()->getName();
                queryLength = topIt2->getSequence()->getSequenceLength();
                queryStart = topIt2->getStartPosition();
                queryEnd = topIt2->getEndPosition();
                targetName = fullNames ? botIt2->getSequence()->getFullName() : botIt2->getSequence()->getName();
                targetLength = botIt2->getSequence()->getSequenceLength();
                targetStart = botIt2->bseg()->getStartPosition();
                targetEnd = botIt2->bseg()->getEndPosition();
                matches = topIt2->getLength();
                snps = countSnps(topIt2, botIt2);
                gaps = 0;
                runningMatch = topIt2->getLength();
            } 
        } else if (found_match) {
            // dump out a running match before we add to the cigar
            if ((cat == 'i' || cat == 'd') && runningMatch > 0) {
                cigar.push_back(std::to_string(runningMatch) + "M");
                runningMatch = 0;
            }
            // extend current paf line
            if (cat == 'i') {
                // extend with inseriton
                hal_index_t ins_len = topIt2->getStartPosition() - topIt1->getEndPosition() - 1;
                assert(ins_len > 0);
                cigar.push_back(std::to_string(ins_len) + "I");
                gaps += ins_len;
            } else if (cat == 'd') {
                // extend with deltion
                hal_index_t del_len;
                if (botIt1->getReversed()) {
                    del_len = botIt1->bseg()->getStartPosition() - botIt2->bseg()->getEndPosition() - 1;
                } else {
                    del_len = botIt2->getStartPosition() - botIt1->getEndPosition() - 1; 
                }
                assert(del_len > 0);
                cigar.push_back(std::to_string(del_len) + "D");
                gaps += del_len;
            }
            if (cat != 'o') {
                // no softclips, so always bookended by a match
                matches += topIt2->getLength();
                snps += countSnps(topIt2, botIt2);
                runningMatch += topIt2->getLength();
            }

            // update the range
            queryEnd = topIt2->getEndPosition();
            targetStart = std::min(targetStart, botIt1->bseg()->getStartPosition());
            targetStart = std::min(targetStart, botIt2->bseg()->getStartPosition());
            targetEnd = std::max(targetEnd, botIt1->bseg()->getEndPosition());
            targetEnd = std::max(targetEnd, botIt2->bseg()->getEndPosition());
        }
            
        // scan the first iterators forward
        topIt1->copy(topIt2);
        botIt1->copy(botIt2);
    }
}
