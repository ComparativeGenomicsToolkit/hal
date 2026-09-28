/*
 * HAL format 3 API: errors.
 *
 * Every failure is a std::runtime_error.  The API's own checks (an unknown genome, a
 * position outside a sequence) throw Error, which says what kind of failure it is;
 * failures reading the files may still arrive as a plain std::runtime_error for now.
 * Lookups that can miss return nullptr or -1 instead of throwing.
 */
#ifndef _HAL3_API_ERROR_H
#define _HAL3_API_ERROR_H

#include <stdexcept>
#include <string>

namespace hal3 {
namespace api {

enum class ErrorKind { Io, Format, Checksum, NotFound, Argument, Unsupported };

class Error : public std::runtime_error {
  public:
    Error(ErrorKind kind, const std::string &what) : std::runtime_error(what), _kind(kind) {}
    ErrorKind kind() const { return _kind; }

  private:
    ErrorKind _kind;
};
}
}

#endif
