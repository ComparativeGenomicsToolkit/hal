/*
 * HAL format 3 API: genomes' DNA.
 *
 * A DnaReader reads bases of any genome, from its .2bit, matching the genome's
 * sequences to the .2bit's by name itself.  It keeps a handle per genome, and each
 * handle keeps what it has read of the .2bit's block lists; recycle(), between units of
 * work, reopens those that have served a lot, which bounds that memory.  One DnaReader
 * per thread.
 */
#ifndef _HAL3_API_DNA_H
#define _HAL3_API_DNA_H

#include "hal3Reader.h"

namespace hal3 {
namespace api {

class DnaReader {
  public:
    explicit DnaReader(std::shared_ptr<const Alignment> alignment);
    ~DnaReader();
    DnaReader(const DnaReader &) = delete;
    DnaReader &operator=(const DnaReader &) = delete;

    /* out = bases [start, end) of sequence seq of g, forward strand: A C G T N, lower
     * case where soft-masked unless noMask.  Throws Error(NotFound) if g has no DNA
     * (see Genome::hasDna), Error(Argument) for a range outside the sequence. */
    void fetch(const Genome *g, uint32_t seq, int64_t start, int64_t end, std::string &out, bool noMask = false);
    void recycle();

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

/* The reverse complement of s, in place (other IUPAC codes kept as they are). */
void reverseComplement(std::string &s);
}
}

#endif
