#include "hal.h"

using namespace std;
using namespace hal;

// Mark that all nodes above this one (but not this one) need to be
// updated.
void markAncestorsForUpdate(Alignment *alignment, string node) {
    Genome *parent = NULL;
    try {
        parent = alignment->openGenome(alignment->getParentName(node));
    } catch (const GenomeNotFoundException& ex) {
        return;  // ignore, not sure why, original code did this.
    }
    MetaData *metadata = parent->getMetaData();
    metadata->set("needsUpdate", "true");
    markAncestorsForUpdate(alignment, parent->getName());
    alignment->closeGenome(parent);
}
