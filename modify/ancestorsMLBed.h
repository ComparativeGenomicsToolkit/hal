#include "halBedScanner.h"
extern "C" {
#include "tree_model.h"
}
class AncestorsMLBed : public hal::BedScanner
{
public:
AncestorsMLBed(TreeModel *mod, hal::AlignmentPtr alignment, hal::Genome *genome, std::map<std::string, int> &nameToId, double threshold, bool writeHal, bool printWrites) : _mod(mod), _alignment(alignment), _genome(genome), _nameToId(nameToId), _threshold(threshold), _writeHal(writeHal), _printWrites(printWrites) {};
  void visitLine();
  TreeModel *_mod;
  hal::AlignmentPtr _alignment;
  hal::Genome *_genome;
  std::map<std::string, int> &_nameToId;
  double _threshold;
  bool _writeHal;
  bool _printWrites;
};
