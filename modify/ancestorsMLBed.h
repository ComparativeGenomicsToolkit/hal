#include "halBedScanner.h"
extern "C" {
#include "tree_model.h"
}
// PHAST code defines min, max macros which conflict with the reserved C++ names.
#undef min
#undef max

class AncestorsMLBed : public hal::BedScanner
{
public:
AncestorsMLBed(TreeModel *mod, hal::AlignmentPtr alignment, const hal::Genome *genome, std::map<std::string, int> &nameToId, double threshold, bool printWrites, bool outputPosts) : _mod(mod), _alignment(alignment), _genome(genome), _nameToId(nameToId), _threshold(threshold), _printWrites(printWrites), _outputPosts(outputPosts) {};
  void visitLine();
  TreeModel *_mod;
  hal::AlignmentConstPtr _alignment;
  const hal::Genome *_genome;
  std::map<std::string, int> &_nameToId;
  double _threshold;
  bool _printWrites;
  bool _outputPosts;
};
