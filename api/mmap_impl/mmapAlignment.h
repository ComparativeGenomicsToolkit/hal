#ifndef _MMAPALIGNMENT_H
#define _MMAPALIGNMENT_H
#include "halAlignment.h"

namespace hal {

class MMapAlignment : public Alignment {
    public:
    // Allocate new array and return the offset.
    hal_index_t allocateNewArray(size_t size);
};
}
#endif
