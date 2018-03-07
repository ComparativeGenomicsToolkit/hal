#include "hal.h"

using namespace std;
using namespace hal;

int main(int argc, char** argv)
{
    CLParserPtr optionsParser = hdf5CLParserInstance();
    optionsParser->setDescription("Calculate coverage by sampling bases.");
    optionsParser->addArgument("halFile", "path to hal file to analyze");
    optionsParser->addArgument("refGenome", "genome to calculate coverage on");
    optionsParser->addOption("numSamples",
                             "Number of bases to sample when calculating coverage",
                             1000000);
    optionsParser->addOption("seed",
                             "Random seed (integer)",
                             0);
    string path;
    string refGenome;
    hal_size_t numSamples;
    int64_t seed;
    try
    {
        optionsParser->parseOptions(argc, argv);
        path = optionsParser->getArgument<string>("halFile");
        refGenome = optionsParser->getArgument<string>("refGenome");
        numSamples = optionsParser->getOption<hal_size_t>("numSamples");
        seed = optionsParser->getOption<int64_t>("seed");
    }
    catch(exception& e)
    {
        cerr << e.what() << endl;
        optionsParser->printUsage(cerr);
        exit(1);
    }

    if (seed == 0) {
        // Default seed. Generate a "random" seed based on the time.
        time_t curTime = time(NULL);
        seed = curTime;
    }
    st_randomSeed(seed);

    AlignmentConstPtr alignment = openHalAlignmentReadOnly(path, optionsParser);
    const Genome *ref = alignment->openGenome(refGenome);
    vector<const Genome *> leafGenomes = getLeafGenomes(alignment);

    map<const Genome *, vector<hal_size_t> > coverage;
    for (size_t i = 0; i < leafGenomes.size(); i++) {
        coverage.insert(make_pair(leafGenomes[i], vector<hal_size_t>()));
    }

    hal_size_t maxDepth = 0;

    for (hal_size_t i = 0; i < numSamples; i++) {
        // Sample (with replacement) a random position in the reference genome.
        hal_index_t pos = st_randomInt64(0, ref->getSequenceLength());
        SegmentIteratorConstPtr refSeg = ref->getTopSegmentIterator();
        refSeg->toSite(pos, true);
        assert(refSeg->getLength() == 1);
        for (size_t j = 0; j < leafGenomes.size(); j++) {
            const Genome *leafGenome = leafGenomes[j];
            set<MappedSegmentConstPtr> segments;
            refSeg->getMappedSegments(segments, leafGenome, NULL,
                                      true, 0, NULL, NULL);
            vector<hal_size_t> &histogram = coverage[leafGenome];
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
        }
    }

    cout << "Genome";
    for (hal_size_t i = 0; i < maxDepth; i++) {
        cout << ", sitesCovered" << i + 1 << "Times";
    }
    cout << endl;

    for (map<const Genome *, vector<hal_size_t> >::iterator it = coverage.begin();
         it != coverage.end(); it++) {
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
