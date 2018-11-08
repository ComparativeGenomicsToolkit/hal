#include "halBedScanner.h"
extern "C" {
#include "tree_model.h"
}
// PHAST code defines min, max macros which conflict with the reserved C++ names.
#undef min
#undef max

using namespace hal;

class AncestorsMLBed : public hal::BedScanner
{
public:
    AncestorsMLBed(TreeModel *mod, const Alignment* alignment, const Genome *genome, std::map<std::string, int> &nameToId, double threshold, bool printWrites, bool outputPosts) : _mod(mod), _alignment(alignment), _genome(genome), _nameToId(nameToId), _threshold(threshold), _printWrites(printWrites), _outputPosts(outputPosts) {};
  void visitLine();
  TreeModel *_mod;
  AlignmentConstPtr _alignment;
  const Genome *_genome;
  std::map<std::string, int> &_nameToId;
  double _threshold;
  bool _printWrites;
  bool _outputPosts;
};
