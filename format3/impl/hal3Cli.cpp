/*
 * HAL format 3: command-line helpers.  See hal3Cli.h.
 */
#include "hal3Cli.h"

#include <cerrno>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace hal3 {

int64_t parseInt(const char *option, const char *value) {
    char *end = nullptr;
    errno = 0;
    long long v = strtoll(value, &end, 10);
    if (!*value || *end || errno) throw CliError(string("--") + option + " needs a whole number, not '" + value + "'");
    return v;
}

double parseDouble(const char *option, const char *value) {
    char *end = nullptr;
    errno = 0;
    double v = strtod(value, &end);
    if (!*value || *end || errno) throw CliError(string("--") + option + " needs a number, not '" + value + "'");
    return v;
}

int runCommand(const char *command, const function<int()> &body) {
    try {
        return body();
    } catch (CliError &e) {
        cerr << "hal " << command << ": " << e.what() << " (hal " << command << " --help for usage)" << endl;
    } catch (exception &e) {
        cerr << "hal " << command << ": " << e.what() << endl;
    }
    return 1;
}
}
