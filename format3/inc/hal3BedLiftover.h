/*
 * HAL format 3: liftover of BED files, written as halLiftover writes them.
 *
 * BedLiftover writes BED exactly as halLiftover does, reproducing its BlockLiftover
 * grouping of mapped pieces into output lines (which depends only on coordinates, not
 * on how HAL happened to segment the genomes) and its handling of BED lines of 3 to
 * 12 columns.  PSL output and --coalescenceLimit are not implemented yet.  It reads
 * BED with HAL's own scanner (halBedScanner.h), so it is kept out of the public API.
 */
#ifndef _HAL3BEDLIFTOVER_H
#define _HAL3BEDLIFTOVER_H

#include "hal3Liftover.h"
#include "halBedScanner.h"

#include <iostream>
#include <list>
#include <set>
#include <string>

namespace hal3 {

class BedLiftover : public hal::BedScanner {
  public:
    BedLiftover(const Genome *src, const Genome *tgt, bool dupes, std::ostream *out);
    /* bedType 0: infer from the first line, as halLiftover */
    void convert(std::istream *in, int bedType);

  protected:
    void visitLine() override;

  private:
    typedef std::list<hal::BedLine> BedList;
    const Genome *_src, *_tgt;
    bool _dupes;
    std::ostream *_out;
    int64_t _srcSeq = -1;
    std::set<std::string> _missed;
    BedList _mappedBlocks, _outBedLines;
    CursorCache _cursors;

    void liftBedInterval(BedList &mapped);
    void liftBlockIntervals();
    void assignBlocksToIntervals();
    bool compatible(const hal::BedLine &tgtBed, const hal::BedLine &newBlock) const;
    void flipBlocks(BedList &bedList) const;
    void cleanResults();
};
}

#endif
