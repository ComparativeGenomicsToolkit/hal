#include "hal.h"
#include "halCLParser.h"

using namespace std;
using namespace hal;

int main(int argc, char **argv) {
    CLParser optionsParser;
    optionsParser.setDescription("Calculate % identity by sampling bases.");
    optionsParser.addArgument("halFile", "path to hal file to analyze");
    optionsParser.addArgument("refGenome", "genome to calculate coverage on");
    optionsParser.addOption("numSamples", "Number of bases to sample when calculating % ID", 1000000);
    optionsParser.addOption("seed", "Random seed (integer)", 0);
    string path;
    string refGenome;
    hal_size_t numSamples;
    int64_t seed;
    try {
        optionsParser.parseOptions(argc, argv);
        path = optionsParser.getArgument<string>("halFile");
        refGenome = optionsParser.getArgument<string>("refGenome");
        numSamples = optionsParser.getOption<hal_size_t>("numSamples");
        seed = optionsParser.getOption<int64_t>("seed");
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

    // Genome -> <# of identical bases to ref, # of total bases aligned to ref>
    map<const Genome *, pair<hal_size_t, hal_size_t>> idStats;
    for (size_t i = 0; i < leafGenomes.size(); i++) {
        idStats.insert(make_pair(leafGenomes[i], make_pair(0, 0)));
    }

    for (hal_size_t i = 0; i < numSamples; i++) {
        // Sample (with replacement) a random position in the reference genome.
        hal_index_t pos = st_randomInt64(0, ref->getSequenceLength());
        SegmentIteratorPtr refSeg = ref->getTopSegmentIterator();
        refSeg->toSite(pos, true);
        assert(refSeg->getLength() == 1);
        string refString;
        refSeg->getString(refString);
        assert(refString.size() == 1);
        if (toupper(refString[0]) == 'N') {
            continue;
        }
        for (size_t j = 0; j < leafGenomes.size(); j++) {
            const Genome *leafGenome = leafGenomes[j];
            MappedSegmentSet segments;
            halMapSegmentSP(refSeg, segments, leafGenome, NULL, true, 0, NULL, NULL);
            if (segments.size() == 1) {
                auto i = segments.begin();
                string tgtString;
                (*i)->getString(tgtString);
                if (toupper(tgtString[0]) == 'N') {
                    continue;
                }
                bool identical = toupper(refString[0]) == toupper(tgtString[0]);
                auto mapIt = idStats.find(leafGenome);
                if (identical) {
                    mapIt->second.first++;
                }
                mapIt->second.second++;
            }
        }
    }

    cout << "Genome, IdenticalSites, AlignedSites, PercentIdentity" << endl;

    for (map<const Genome *, pair<hal_size_t, hal_size_t>>::iterator it = idStats.begin(); it != idStats.end(); it++) {
        string name = it->first->getName();
        hal_size_t identicalBases = it->second.first;
        hal_size_t alignedBases = it->second.second;
        cout << name << ", "
             << identicalBases << ", "
             << alignedBases << ", "
             << 100.0 * ((double) identicalBases) / alignedBases
             << endl;
    }
}
