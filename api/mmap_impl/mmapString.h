#ifndef _MMAP_STRING_H
#define _MMAP_STRING_H
#include "mmapArray.h"
namespace hal {
class MMapString : public MMapArray<char> {
    public:
    MMapString(MMapAlignment *alignment, const std::string &string) : MMapArray(alignment), _string(string) { set(_string); };
    MMapString(MMapAlignment *alignment, size_t offset) : MMapArray(alignment, offset) { read(); };
    const char *c_str() { return getSlice(0, getLength() - 1); }
    const std::string &get() { return _string; }
    size_t set(const std::string &string) {
        _string = string;
        setLength(string.size() + 1);
        for (size_t i = 0; i < string.size(); i++) {
            *(this->operator[](i)) = string[i];
        }
        *(this->operator[](string.size())) = '\0';
        return getOffset();
    }
    private:
    void read() {
        _string = c_str();
    };
    std::string _string;
};
}
#endif
