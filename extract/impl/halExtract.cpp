/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hal.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace hal;

static void getDimensions(AlignmentConstPtr outAlignment, const Genome *genome, vector<Sequence::Info> &dimensions);

static void copyGenome(const Genome *inGenome, Genome *outGenome);

static void extractTree(AlignmentConstPtr inAlignment, AlignmentPtr outAlignment, const string &rootName);

static void extract(AlignmentConstPtr inAlignment, AlignmentPtr outAlignment, const string &rootName);

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("inHalPath", "input hal file");
    optionsParser.addArgument("outHalPath", "output hal file");
    optionsParser.addOption("outputFormat", "format for output hal file (same as input file by default)", "");
    optionsParser.addOption("root", "root of subtree to extract", "\"\"");
}

int main(int argc, char **argv) {
    CLParser optionsParser(CREATE_ACCESS);
    initParser(optionsParser);

    string inHalPath;
    string outHalPath;
    string rootName;
    string outputFormat;
    try {
        optionsParser.parseOptions(argc, argv);
        inHalPath = optionsParser.getArgument<string>("inHalPath");
        outHalPath = optionsParser.getArgument<string>("outHalPath");
        rootName = optionsParser.getOption<string>("root");
        outputFormat = optionsParser.getOption<string>("outputFormat");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    try {
        AlignmentConstPtr inAlignment(openHalAlignment(inHalPath, &optionsParser));
        if (inAlignment->getNumGenomes() == 0) {
            throw hal_exception("input hal alignmenet is empty");
        }

        if (outputFormat.empty()) {
            // No alignment format specified, just use the same as the input format.
            outputFormat = inAlignment->getStorageFormat();
        }

        AlignmentPtr outAlignment(
            openHalAlignment(outHalPath, &optionsParser, READ_ACCESS | WRITE_ACCESS | CREATE_ACCESS, outputFormat));
        if (outAlignment->getNumGenomes() != 0) {
            throw hal_exception("output hal alignmenet cannot be initialized");
        }

        if (rootName == "\"\"" || inAlignment->getNumGenomes() == 0) {
            rootName = inAlignment->getRootName();
        }

        extractTree(inAlignment, outAlignment, rootName);
        extract(inAlignment, outAlignment, rootName);
        outAlignment->close();
    } catch (hal_exception &e) {
        cerr << "hal exception caught: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
    return 0;
}

void getDimensions(AlignmentConstPtr outAlignment, const Genome *genome, vector<Sequence::Info> &dimensions) {
    assert(dimensions.size() == 0);

    bool root = outAlignment->getParentName(genome->getName()).empty();
    bool leaf = outAlignment->getChildNames(genome->getName()).empty();

    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        Sequence::Info info(sequence->getName(), sequence->getSequenceLength(), root ? 0 : sequence->getNumTopSegments(),
                            leaf ? 0 : sequence->getNumBottomSegments());
        dimensions.push_back(info);
    }
}

void copyGenome(const Genome *inGenome, Genome *outGenome) {
    DnaIteratorPtr inDna = inGenome->getDnaIterator();
    DnaIteratorPtr outDna = outGenome->getDnaIterator();
    hal_size_t n = inGenome->getSequenceLength();
    assert(n == outGenome->getSequenceLength());
    for (; (hal_size_t)inDna->getArrayIndex() < n; inDna->toRight(), outDna->toRight()) {
        outDna->setBase(inDna->getBase());
    }
    outDna->flush();

    TopSegmentIteratorPtr inTop = inGenome->getTopSegmentIterator();
    TopSegmentIteratorPtr outTop = outGenome->getTopSegmentIterator();
    n = outGenome->getNumTopSegments();
    assert(n == 0 || n == inGenome->getNumTopSegments());
    for (; (hal_size_t)inTop->getArrayIndex() < n; inTop->toRight(), outTop->toRight()) {
        outTop->setCoordinates(inTop->getStartPosition(), inTop->getLength());
        outTop->tseg()->setParentIndex(inTop->tseg()->getParentIndex());
        outTop->tseg()->setParentReversed(inTop->tseg()->getParentReversed());
        outTop->tseg()->setBottomParseIndex(inTop->tseg()->getBottomParseIndex());
        outTop->tseg()->setNextParalogyIndex(inTop->tseg()->getNextParalogyIndex());
    }

    BottomSegmentIteratorPtr inBot = inGenome->getBottomSegmentIterator();
    BottomSegmentIteratorPtr outBot = outGenome->getBottomSegmentIterator();
    n = outGenome->getNumBottomSegments();
    assert(n == 0 || n == inGenome->getNumBottomSegments());
    hal_size_t nc = inGenome->getNumChildren();
    assert(nc == outGenome->getNumChildren());
    for (; (hal_size_t)inBot->getArrayIndex() < n; inBot->toRight(), outBot->toRight()) {
        outBot->setCoordinates(inBot->getStartPosition(), inBot->getLength());
        for (hal_size_t child = 0; child < nc; ++child) {
            outBot->bseg()->setChildIndex(child, inBot->bseg()->getChildIndex(child));
            outBot->bseg()->setChildReversed(child, inBot->bseg()->getChildReversed(child));
        }
        if (outGenome->getAlignment()->getRootName() == outGenome->getName()) {
            outBot->bseg()->setTopParseIndex(NULL_INDEX);
        } else {
            outBot->bseg()->setTopParseIndex(inBot->bseg()->getTopParseIndex());
        }
    }

    const map<string, string> &meta = inGenome->getMetaData()->getMap();
    map<string, string>::const_iterator i = meta.begin();
    for (; i != meta.end(); ++i) {
        outGenome->getMetaData()->set(i->first, i->second);
    }
}

static void extractTree(AlignmentConstPtr inAlignment, AlignmentPtr outAlignment, const string &rootName) {
    const Genome *genome = inAlignment->openGenome(rootName);
    if (genome == NULL) {
        throw hal_exception(string("Genome not found: ") + rootName);
    }
    Genome *newGenome = NULL;
    if (outAlignment->getNumGenomes() == 0 || genome->getParent() == NULL) {
        newGenome = outAlignment->addRootGenome(rootName);
    } else {
        const Genome *parent = genome->getParent();
        assert(parent != NULL);
        newGenome =
            outAlignment->addLeafGenome(rootName, parent->getName(), inAlignment->getBranchLength(parent->getName(), rootName));
    }
    assert(newGenome != NULL);

    inAlignment->closeGenome(genome);
    outAlignment->closeGenome(newGenome);

    vector<string> childNames = inAlignment->getChildNames(rootName);
    for (size_t i = 0; i < childNames.size(); ++i) {
        extractTree(inAlignment, outAlignment, childNames[i]);
    }
}

void extract(AlignmentConstPtr inAlignment, AlignmentPtr outAlignment, const string &rootName) {
    const Genome *genome = inAlignment->openGenome(rootName);
    Genome *newGenome = outAlignment->openGenome(rootName);
    assert(newGenome != NULL);

    vector<Sequence::Info> dimensions;
    getDimensions(inAlignment, genome, dimensions);
    newGenome->setDimensions(dimensions);

    cout << "Extracting " << genome->getName() << endl;
    copyGenome(genome, newGenome);

    inAlignment->closeGenome(genome);
    outAlignment->closeGenome(newGenome);

    vector<string> childNames = inAlignment->getChildNames(rootName);
    for (size_t i = 0; i < childNames.size(); ++i) {
        extract(inAlignment, outAlignment, childNames[i]);
    }
}
