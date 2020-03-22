/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "hal.h"
#include "halCLParser.h"

#include "hal2psl.h"
#include "psl_io.h"
#include "psl_merger.h"

using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("alignment", "input file in HAL or PSL format (PSL must specify --alignmentIsPsl)");
    optionsParser.addOptionFlag("alignmentIsPsl", "alignment is in PSL format", false);
    optionsParser.addOption("queryGenome", "source genome", "\"\"");
    optionsParser.addOption("targetGenome", "reference genome name", "\"\"");
    optionsParser.addArgument("outPslPath", "output psl file ffor synteny blocks");
    optionsParser.addOption("minBlockSize", "lower bound on synteny block length", 5000);
    optionsParser.addOption("maxAnchorDistance", "upper bound on distance for syntenic psl blocks", 5000);
    optionsParser.addOption("queryChromosome", "chromosome to infer synteny (default is whole genome)", "\"\"");
    optionsParser.setDescription("Convert alignments into synteny blocks");
}

const Genome *openGenomeOrThrow(const Alignment *alignment, const std::string &genomeName) {
    const Genome *genome = alignment->openGenome(genomeName);
    if (genome == NULL) {
        throw hal_exception(std::string("Reference genome, ") + genomeName + ", not found in alignment");
    }
    return genome;
}

const Alignment *openAlignmentOrThrow(const std::string &alignmentFile, const CLParser &optionsParser) {

    auto alignment = openHalAlignment(alignmentFile, &optionsParser);
    if (alignment->getNumGenomes() == 0) {
        throw hal_exception("hal alignment is empty");
    }
    return alignment;
}

void validateInputOrThrow(const std::string &queryGenomeName, const std::string &targetGenomeName,
                          const std::string &queryChromosome) {
    if (queryGenomeName == "\"\"" || targetGenomeName == "\"\"") {
        throw hal_exception("--queryGenome and --targetGenome and --queryChromosome must be"
                            "specified");
    }

    if (queryGenomeName == targetGenomeName) {
        throw hal_exception("--queryGenome and --targetGenome must be"
                            "different");
    }
}

static void makeSyntenyBlocks(std::vector<PslBlock>& blocks, hal_size_t minBlockSize,
                              hal_size_t maxAnchorDistance, std::ofstream &pslFh) {
    auto merged_blocks = dag_merge(blocks, minBlockSize, maxAnchorDistance);
    psl_io::write_psl(merged_blocks, pslFh);
}

static void syntenyFromPsl(std::string alignmentFile, hal_size_t minBlockSize,
                           hal_size_t maxAnchorDistance, std::string outPslPath) {
    auto blocks = psl_io::get_blocks_set(alignmentFile);
    std::ofstream pslFh;
    pslFh.exceptions(std::ofstream::failbit|std::ofstream::badbit);
    pslFh.open(outPslPath, std::ofstream::out);
    makeSyntenyBlocks(blocks, minBlockSize, maxAnchorDistance, pslFh);
    pslFh.close();
}

static std::vector<std::string> getChromNames(const Genome *genome) {
    std::vector<std::string> chromNames;
    for (auto seqIt = genome->getSequenceIterator();
         not seqIt->atEnd(); seqIt->toNext()) {
        chromNames.push_back(seqIt->getSequence()->getName());
    }
    std::sort(chromNames.begin(), chromNames.end());
    return chromNames;
}

static void syntenyBlockForChrom(const Alignment *alignment,
                                 const Genome *targetGenome, const Genome *queryGenome,
                                 std::string queryChromosome, hal_size_t minBlockSize,
                                 hal_size_t maxAnchorDistance, std::ofstream &pslFh) {
    auto hal2psl = hal::Hal2Psl();
    auto blocks = hal2psl.convert2psl(alignment, queryGenome, targetGenome, queryChromosome);
    makeSyntenyBlocks(blocks, minBlockSize, maxAnchorDistance, pslFh);
}


/* do one chromosome at a time to reduce memory */
static void syntenyFromHal(const Alignment *alignment, std::string queryGenomeName,
                           std::string targetGenomeName, std::string queryChromosome,
                           hal_size_t minBlockSize, hal_size_t maxAnchorDistance, std::string outPslPath) {
    auto targetGenome = openGenomeOrThrow(alignment, targetGenomeName);
    auto queryGenome = openGenomeOrThrow(alignment, queryGenomeName);
    std::vector<std::string> chromNames;
    if (queryChromosome != "\"\"") {
        chromNames.push_back(queryChromosome);
    } else {
        chromNames = getChromNames(queryGenome);
    }


    std::ofstream pslFh;
    pslFh.exceptions(std::ofstream::failbit|std::ofstream::badbit);
    pslFh.open(outPslPath, std::ofstream::out);
    for (auto chromIt = chromNames.begin(); chromIt != chromNames.end(); chromIt++) {
        syntenyBlockForChrom(alignment, targetGenome, queryGenome,
                             *chromIt, minBlockSize, maxAnchorDistance, pslFh);
    }
    pslFh.close();
}

int main(int argc, char *argv[]) {
    CLParser optionsParser;
    initParser(optionsParser);
    std::string alignmentFile;
    bool alignmentIsPsl;
    std::string outPslPath;
    std::string queryGenomeName;
    std::string targetGenomeName;
    std::string queryChromosome;
    hal_size_t minBlockSize;
    hal_size_t maxAnchorDistance;
    try {
        optionsParser.parseOptions(argc, argv);
        alignmentFile = optionsParser.getArgument<std::string>("alignment");
        alignmentIsPsl = optionsParser.getFlag("alignmentIsPsl");
        queryGenomeName = optionsParser.getOption<std::string>("queryGenome");
        targetGenomeName = optionsParser.getOption<std::string>("targetGenome");
        outPslPath = optionsParser.getArgument<std::string>("outPslPath");
        minBlockSize = optionsParser.getOption<hal_size_t>("minBlockSize");
        maxAnchorDistance = optionsParser.getOption<hal_size_t>("maxAnchorDistance");
        queryChromosome = optionsParser.getOption<std::string>("queryChromosome");
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        optionsParser.printUsage(std::cerr);
        exit(1);
    }

    validateInputOrThrow(queryGenomeName, targetGenomeName, queryChromosome);
    try {
        std::vector<PslBlock> blocks;
        if (alignmentIsPsl) {
            syntenyFromPsl(alignmentFile, minBlockSize, maxAnchorDistance, outPslPath);
        } else {
            auto alignment = openAlignmentOrThrow(alignmentFile, optionsParser);
            syntenyFromHal(alignment, queryGenomeName, targetGenomeName, queryChromosome, minBlockSize, maxAnchorDistance, outPslPath);
            alignment->close();
        }
    } catch (std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
