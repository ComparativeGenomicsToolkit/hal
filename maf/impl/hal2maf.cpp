/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafBed.h"
#include "halMafExport.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("halFile", "input hal file");
    optionsParser.addArgument("mafFile", "output maf file (or \"stdout\" to "
                                         "pipe to standard output)");
    optionsParser.addOption("refGenome", "name of reference genome (root if empty)", "");
    optionsParser.addOption("refSequence", "name of reference sequence within reference genome"
                                           " (all sequences if empty)",
                            "");
    optionsParser.addOption("refTargets", "bed file coordinates of intervals in the reference "
                                          "genome to export (or \"stdin\" to pipe from "
                                          "standard input)",
                            "");
    optionsParser.addOption("start", "coordinate within the sequence, requires --refSequence", 0);
    optionsParser.addOption("length", "length of the sequence, requires --refSequence. "
                                      "If set to 0, the entire genome is converted",
                            0);
    optionsParser.addOption("rootGenome", "name of root genome (none if empty)", "");
    optionsParser.addOption("targetGenomes", "comma-separated (no spaces) list of target genomes "
                                             "(others are excluded) (vist all if empty)",
                            "");
    optionsParser.addOption("maxRefGap", "maximum gap length in reference", 0);
    optionsParser.addOptionFlag("noDupes", "ignore paralogy edges", false);
    optionsParser.addOptionFlag("noAncestors", "don't write ancestral sequences. IMPORTANT: "
                                               "Must be used in conjunction with --refGenome"
                                               " to set a non-ancestral genome as the reference"
                                               " because the default reference is the root.",
                                false);
    optionsParser.addOptionFlag("onlySequenceNames", "use only sequence names "
                                                     "for output names.  By default, the UCSC convention of Genome.Sequence "
                                                     "is used",
                                false);
    optionsParser.addOptionFlag("unique", "only write column whose left-most reference "
                                          "coordinate is in the specified range.  this "
                                          "is used to insure that the same column isnt "
                                          "sampled twice (due to ducplications) by mafs "
                                          "generated on distinct ranges.",
                                false);
    optionsParser.addOptionFlag("append", "append to instead of overwrite output file.", false);
    optionsParser.addOption("maxBlockLen", "maximum length of MAF block in output", MafBlock::defaultMaxLength);
    optionsParser.addOptionFlag("global", "output all columns in alignment, "
                                          "ignoring refGenome, refSequence, etc. flags",
                                false);
    optionsParser.addOptionFlag("printTree", "print a gene tree for every block", false);
    optionsParser.addOptionFlag("onlyOrthologs", "make only orthologs to the "
                                                 "reference appear in the MAF blocks",
                                false);

    optionsParser.setDescription("Convert hal database to maf.");
}

/* Parsed options */
class MafOptions {
  public:
    string halPath;
    string mafPath;
    string refGenomeName;
    string rootGenomeName;
    string targetGenomes;
    string refSequenceName;
    string refTargetsPath;
    hal_index_t start;
    hal_size_t length;
    hal_size_t maxRefGap;
    bool noDupes;
    bool noAncestors;
    bool ucscNames;
    bool unique;
    bool append;
    bool global;
    bool printTree;
    bool onlyOrthologs;
    hal_index_t maxBlockLen;
};

/* This empty string options specified using the old convention of '""' rather than
 * just an empty string. FIXME: this should be removed. */
static string fixString(const string &s) {
    if (s == "\"\"") {
        cerr << "WARNING missing string arguments should be specified as empty strings, not the obsolete '\"\"'" << endl;
        return "";
    } else {
        return s;
    }
}

static void hal2mafWithTargets(const MafOptions &opts, AlignmentConstPtr alignment, const Genome *refGenome,
                               set<const Genome *> &targetSet, MafExport &mafExport, ostream &mafStream) {
    ifstream bedFileStream;
    ifstream refTargetsStream;
    if (opts.refTargetsPath != "stdin") {
        bedFileStream.open(opts.refTargetsPath);
        if (!refTargetsStream) {
            throw hal_exception("Error opening " + opts.refTargetsPath);
        }
    }
    istream &bedStream = opts.refTargetsPath != "stdin" ? bedFileStream : cin;
    MafBed mafBed(mafStream, alignment, refGenome, targetSet, mafExport);
    mafBed.scan(&bedStream);
}

static void hal2maf(AlignmentConstPtr alignment, const MafOptions &opts) {
    const Genome *rootGenome = NULL;
    set<const Genome *> targetSet;
    if (opts.rootGenomeName != "") {
        rootGenome = alignment->openGenome(opts.rootGenomeName);
        if (rootGenome == NULL) {
            throw hal_exception("Root genome " + opts.rootGenomeName + ", not found in alignment");
        }
        if (opts.rootGenomeName != alignment->getRootName()) {
            getGenomesInSubTree(rootGenome, targetSet);
        }
    }

    if (opts.targetGenomes != "") {
        vector<string> targetNames = chopString(opts.targetGenomes, ",");
        for (size_t i = 0; i < targetNames.size(); ++i) {
            const Genome *tgtGenome = alignment->openGenome(targetNames[i]);
            if (tgtGenome == NULL) {
                throw hal_exception(string("Target genome, ") + targetNames[i] + ", not found in alignment");
            }
            targetSet.insert(tgtGenome);
        }
    }

    const Genome *refGenome = NULL;
    if (opts.refGenomeName != "") {
        refGenome = alignment->openGenome(opts.refGenomeName);
        if (refGenome == NULL) {
            throw hal_exception("Reference genome, " + opts.refGenomeName + ", not found in alignment");
        }
    } else {
        refGenome = alignment->openGenome(alignment->getRootName());
    }
    if (opts.noAncestors && refGenome->getNumChildren() != 0 && !opts.global) {
        throw hal_exception(string("Since the reference genome to be used for the"
                                   " MAF is ancestral (") +
                            refGenome->getName() + "), the --noAncestors option is invalid.  The "
                                                   "--refGenome option can be used to specify a "
                                                   "different reference.");
    }

    const Sequence *refSequence = NULL;
    if (opts.refSequenceName != "") {
        refSequence = refGenome->getSequence(opts.refSequenceName);
        if (refSequence == NULL) {
            throw hal_exception(string("Reference sequence, ") + opts.refSequenceName + ", not found in reference genome, " +
                                refGenome->getName());
        }
    }

    ios_base::openmode openFlags = ios_base::out;
    if (opts.append == true) {
        openFlags |= ios_base::app;
    }
    ofstream mafFileStream;
    if (opts.mafPath != "stdout") {
        mafFileStream.open(opts.mafPath, openFlags);
        if (!mafFileStream) {
            throw hal_exception("Error opening " + opts.mafPath);
        }
    }
    ostream &mafStream = opts.mafPath != "stdout" ? mafFileStream : cout;

    MafExport mafExport;
    mafExport.setMaxRefGap(opts.maxRefGap);
    mafExport.setNoDupes(opts.noDupes);
    mafExport.setNoAncestors(opts.noAncestors);
    mafExport.setUcscNames(opts.ucscNames);
    mafExport.setUnique(opts.unique);
    mafExport.setAppend(opts.append);
    mafExport.setMaxBlockLength(opts.maxBlockLen);
    mafExport.setPrintTree(opts.printTree);
    mafExport.setOnlyOrthologs(opts.onlyOrthologs);

    if (opts.refTargetsPath != "") {
        hal2mafWithTargets(opts, alignment, refGenome, targetSet, mafExport, mafStream);
    } else if (opts.global) {
        mafExport.convertEntireAlignment(mafStream, alignment);
    } else if (refSequence != NULL) {
        mafExport.convertSequence(mafStream, alignment, refSequence, opts.start, opts.length, targetSet);
    } else {
        for (SequenceIteratorPtr seqIt(refGenome->getSequenceIterator()); not seqIt->atEnd(); seqIt->toNext()) {
            mafExport.convertSequence(mafStream, alignment, seqIt->getSequence(), opts.start, opts.length, targetSet);
        }
    }
    if (opts.mafPath != "stdout") {
        // dont want to leave a size 0 file when there's not ouput because
        // it can make some scripts (ie that process a maf for each contig)
        // obnoxious (presently the case for halPhlyoPTrain which uses
        // hal2mafMP --splitBySequence). FIXME: this can also break stuff that
        // has dependencies, so drop it.
        if (mafFileStream.tellp() == (streampos)0) {
            std::remove(opts.mafPath.c_str());
        }
    }
}

int main(int argc, char **argv) {
    CLParser optionsParser;
    initParser(optionsParser);
    MafOptions opts;
    try {
        optionsParser.parseOptions(argc, argv);
        opts.halPath = fixString(optionsParser.getArgument<string>("halFile"));
        opts.mafPath = fixString(optionsParser.getArgument<string>("mafFile"));
        opts.refGenomeName = fixString(optionsParser.getOption<string>("refGenome"));
        opts.rootGenomeName = fixString(optionsParser.getOption<string>("rootGenome"));
        opts.targetGenomes = fixString(optionsParser.getOption<string>("targetGenomes"));
        opts.refSequenceName = fixString(optionsParser.getOption<string>("refSequence"));
        opts.refTargetsPath = fixString(optionsParser.getOption<string>("refTargets"));
        opts.start = optionsParser.getOption<hal_index_t>("start");
        opts.length = optionsParser.getOption<hal_size_t>("length");
        opts.maxRefGap = optionsParser.getOption<hal_size_t>("maxRefGap");
        opts.noDupes = optionsParser.getFlag("noDupes");
        opts.noAncestors = optionsParser.getFlag("noAncestors");
        opts.ucscNames = !optionsParser.getFlag("onlySequenceNames");
        opts.unique = optionsParser.getFlag("unique");
        opts.append = optionsParser.getFlag("append");
        opts.global = optionsParser.getFlag("global");
        opts.printTree = optionsParser.getFlag("printTree");
        opts.maxBlockLen = optionsParser.getOption<hal_index_t>("maxBlockLen");
        opts.onlyOrthologs = optionsParser.getFlag("onlyOrthologs");

        if (((opts.length != 0) || (opts.start != 0)) && (opts.refSequenceName == "")) {
            throw hal_exception("--start and --length require --refSequenceName");
        }
        if (opts.rootGenomeName != "" && opts.targetGenomes != "") {
            throw hal_exception("--rootGenome and --targetGenomes options are "
                                "mutually exclusive");
        }
        if ((not opts.refTargetsPath.empty()) and
            ((opts.start != 0) || (opts.length != 0) || (not opts.refSequenceName.empty()))) {
            throw hal_exception("--refSequence, --start, and --length options are unsupported when using BED input");
        }

    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
    try {
        AlignmentConstPtr alignment(openHalAlignment(opts.halPath, &optionsParser));
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("hal alignmenet is empty");
        }

        hal2maf(alignment, opts);
    } catch (hal_exception &e) {
        cerr << "hal exception caught: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }

    return 0;
}
