#include "halGlobalColumnIterator.h"
#include "halColumnIterator.h"
#include "hal.h"

using namespace std;
using namespace hal;

GlobalColumnIterator::GlobalColumnIterator(AlignmentConstPtr alignment,
                                           const set<const Genome *> *targets,
                                           bool noDupes, bool noAncestors,
                                           bool reverseStrand) :
    _noDupes(noDupes),
    _noAncestors(noAncestors),
    _reverseStrand(reverseStrand)
{
    if (targets == NULL) {
        // Use all leaf genomes as targets by default.
        vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
        for (hal_size_t i = 0; i < leafNames.size(); i++) {
            const Genome *genome = alignment->openGenome(leafNames[i]);
            assert(genome != NULL);
            _targets.insert(genome);
        }
    } else {
        _targets = *targets;
    }
    _curGenome = _targets.begin();
    _genomeColIt = (*_curGenome)->getColumnIterator(&_targets, 0, 0,
                                                    NULL_INDEX, _noDupes,
                                                    _noAncestors,
                                                    _reverseStrand);
}

GlobalColumnIterator::~GlobalColumnIterator()
{
    clearVisitCache();
}

void GlobalColumnIterator::toRight()
{
    if (_genomeColIt->lastColumn()) {
        // Need to switch to the next genome in the target set.

        // Copy over the updated visit cache information. This is a
        // deep copy, so it's slow, but necessary to preserve the
        // column iterator ownership of the visit cache
        clearVisitCache();
        ColumnIterator::VisitCache *newVisitCache = _genomeColIt->getVisitCache();
        for(ColumnIterator::VisitCache::iterator it = newVisitCache->begin();
            it != newVisitCache->end(); it++) {
            _globalVisitCache[it->first] = new PositionCache(*it->second);
        }

        // We don't have to worry about whether there is another
        // genome in the set because column iterators are allowed to
        // crash if lastColumn() is true and toRight() is called.
        _curGenome++;
        _genomeColIt = (*_curGenome)->getColumnIterator(&_targets, 0, 0,
                                                        NULL_INDEX, _noDupes,
                                                        _noAncestors,
                                                        _reverseStrand);
        _genomeColIt->setVisitCache(&_globalVisitCache);
    } else {
        // Can safely go to the next unvisited site in the genome.
        _genomeColIt->toRight();
    }
}

const ColumnIterator::ColumnMap *GlobalColumnIterator::getColumnMap() const
{
    return _genomeColIt->getColumnMap();
}

bool GlobalColumnIterator::lastColumn() const
{
    set<const Genome *>::const_iterator tmpIt = _curGenome;
    tmpIt++;
    return _genomeColIt->lastColumn() && tmpIt == _targets.end();
}

void GlobalColumnIterator::clearVisitCache()
{
    for (VisitCache::iterator i = _globalVisitCache.begin();
         i != _globalVisitCache.end(); ++i)
    {
        delete i->second;
    }
    _globalVisitCache.clear();
}
