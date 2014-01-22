#include "hal.h"

using namespace std;
using namespace hal;

// Mark that all nodes above this one (but not this one) need to be
// updated.
void markAncestorsForUpdate(AlignmentPtr alignment, string node)
{
  Genome *parent = alignment->openGenome(alignment->getParentName(node));
  if (!parent) {
    return;
  }
  MetaData *metadata = parent->getMetaData();
  metadata->set("needsUpdate", "true");
  markAncestorsForUpdate(alignment, parent->getName());
  alignment->closeGenome(parent);
}
