#include "hal.h"
#include "halCLParser.h"

using namespace std;
using namespace hal;

int main(int argc, char **argv) {
    CLParser optionsParser;
    optionsParser.setDescription("Calculate coverage by sampling bases.");
    optionsParser.addArgument("halFile", "path to hal file to analyze");
    optionsParser.addArgument("refGenome", "genome to calculate coverage on");
    optionsParser.addOption("numSamples", "Number of bases to sample when calculating coverage", 1000000);
    optionsParser.addOption("seed", "Random seed (integer)", 0);
    optionsParser.addOptionFlag("bySequence", "provide coverage breakdown by sequence in reference genome", false);
    
    string path;
    string refGenome;
    hal_size_t numSamples;
    int64_t seed;
    bool bySequence;
    try {
        optionsParser.parseOptions(argc, argv);
        path = optionsParser.getArgument<string>("halFile");
        refGenome = optionsParser.getArgument<string>("refGenome");
        numSamples = optionsParser.getOption<hal_size_t>("numSamples");
        seed = optionsParser.getOption<int64_t>("seed");
        bySequence = optionsParser.getFlag("bySequence");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    if (seed == 0) {
        // Default seed. Generate a "random" seed based on the time.
        time_t curTime = time(NULL);
        seed = curTime;
    }
    st_randomSeed(seed);

    AlignmentConstPtr alignment(openHalAlignment(path, &optionsParser));
    const Genome *ref = alignment->openGenome(refGenome);
    vector<const Genome *> leafGenomes = getLeafGenomes(alignment.get());

    // null key means whole genome here
    map<const Sequence*, map<const Genome*, vector<hal_size_t>>> coverage_by_sequence;
    vector<const Sequence*> sequences = {NULL};
    if (bySequence) {
        for (SequenceIteratorPtr si = ref->getSequenceIterator(); !si->atEnd(); si->toNext()) {
            sequences.push_back(si->getSequence());
        }
    }    
    for (const Sequence* sequence : sequences) {
        map<const Genome*, vector<hal_size_t>>& coverage = coverage_by_sequence[sequence];
        for (size_t i = 0; i < leafGenomes.size(); i++) {
            coverage.insert(make_pair(leafGenomes[i], vector<hal_size_t>()));
        }
    }

    hal_size_t maxDepth = 0;

    map<const Genome*, vector<hal_size_t>>& genome_coverage = coverage_by_sequence[NULL];
    for (hal_size_t i = 0; i < numSamples; i++) {
        // Sample (with replacement) a random position in the reference genome.
        hal_index_t pos = st_randomInt64(0, ref->getSequenceLength());
        SegmentIteratorPtr refSeg = ref->getTopSegmentIterator();
        refSeg->toSite(pos, true);
        assert(refSeg->getLength() == 1);
        map<const Genome*, vector<hal_size_t>>* sequence_coverage = NULL;
        if (bySequence) {
            sequence_coverage = &coverage_by_sequence[refSeg->getSequence()];
        }
        for (size_t j = 0; j < leafGenomes.size(); j++) {
            const Genome *leafGenome = leafGenomes[j];
            MappedSegmentSet segments;
            halMapSegmentSP(refSeg, segments, leafGenome, NULL, true, 0, NULL, NULL);
            vector<hal_size_t> &histogram = genome_coverage[leafGenome];
            hal_size_t depth = segments.size();
            if (depth > maxDepth) {
                maxDepth = depth;
            }
            if (histogram.size() < depth) {
                histogram.resize(depth, 0);
            }
            for (size_t k = 0; k < depth; k++) {
                histogram[k] += 1;
            }
            if (sequence_coverage) {
                vector<hal_size_t> &seq_histogram = sequence_coverage->at(leafGenome);
                if (seq_histogram.size() < depth) {
                    seq_histogram.resize(depth, 0);
                }
                for (size_t k = 0; k < depth; k++) {
                    seq_histogram[k] += 1;
                }
            }
        }
    }

    cout << "Genome";
    for (hal_size_t i = 0; i < maxDepth; i++) {
        cout << ", sitesCovered" << i + 1 << "Times";
    }
    cout << endl;

    for (const Sequence* refseq : sequences) {
        if (refseq != NULL) {
            cerr << "\nCoverage on " << refseq->getName() << endl;
        }
        map<const Genome*, vector<hal_size_t>>& coverage = coverage_by_sequence[refseq];
        for (map<const Genome *, vector<hal_size_t>>::iterator it = coverage.begin(); it != coverage.end(); it++) {
            string name = it->first->getName();
            cout << name;
            vector<hal_size_t> histogram = it->second;
            for (hal_size_t i = 0; i < maxDepth; i++) {
                if (i < histogram.size()) {
                    cout << ", " << histogram[i];
                } else {
                    cout << ", 0";
                }
            }
            cout << endl;
        }
    }
}
