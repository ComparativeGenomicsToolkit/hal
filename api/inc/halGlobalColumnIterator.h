#ifndef _GLOBALCOLUMNITERATOR_H
#define _GLOBALCOLUMNITERATOR_H

#include "halColumnIterator.h"
#include "halCommon.h"

namespace hal {
class GlobalColumnIterator
{
public:
    // FIXME: get rid of these after global/local colIt is done properly
    typedef std::map<const Genome*, PositionCache*> VisitCache;
    GlobalColumnIterator(AlignmentConstPtr alignment,
                         const std::set<const Genome *> *targets, bool noDupes,
                         bool noAncestors, bool reverseStrand);
    ~GlobalColumnIterator();

    virtual void toRight();
    virtual bool lastColumn() const;
    virtual const ColumnIterator::ColumnMap *getColumnMap() const;
//    virtual void defragment();
//    virtual void print(std::ostream &os) const;
private:
    void clearVisitCache();
    ColumnIteratorConstPtr _genomeColIt;
    VisitCache _globalVisitCache;
    std::set<const Genome *> _targets;
    std::set<const Genome *>::const_iterator _curGenome;
    bool _noDupes;
    bool _noAncestors;
    bool _reverseStrand;
};
}
#endif
