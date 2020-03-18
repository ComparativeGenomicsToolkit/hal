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
    optionsParser.addArgument("alignment", "input file in HAL or PSL format (PSL must end in .psl)");
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

int main(int argc, char *argv[]) {
    CLParser optionsParser;
    initParser(optionsParser);
    std::string alignmentFile;
    std::string outPslPath;
    std::string queryGenomeName;
    std::string targetGenomeName;
    std::string queryChromosome;
    hal_size_t minBlockSize;
    hal_size_t maxAnchorDistance;
    try {
        optionsParser.parseOptions(argc, argv);
        alignmentFile = optionsParser.getArgument<std::string>("alignment");
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
        if (alignmentFile.substr(alignmentFile.find_last_of(".") + 1) == "psl") {
            blocks = psl_io::get_blocks_set(alignmentFile);
        } else {
            auto alignment = openAlignmentOrThrow(alignmentFile, optionsParser);
            auto targetGenome = openGenomeOrThrow(alignment, targetGenomeName);
            auto queryGenome = openGenomeOrThrow(alignment, queryGenomeName);

            auto hal2psl = hal::Hal2Psl();
            blocks = hal2psl.convert2psl(alignment, queryGenome, targetGenome, queryChromosome);
        }
        std::cout << "merging " << blocks.size() << " blocks" << std::endl;
        auto merged_blocks = dag_merge(blocks, minBlockSize, maxAnchorDistance);
        std::cout << "writing psl " << std::endl;
        psl_io::write_psl(merged_blocks, outPslPath);
    } catch (std::exception &e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
