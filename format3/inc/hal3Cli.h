/*
 * HAL format 3: the `hal` command line.
 *
 * One binary, `hal <command> [options] <args>`, as vg and samtools.  Each command is a
 * function main_<name>(argc, argv) with argv[0] its name, parsing its options with
 * getopt_long; positional arguments keep the order of the HAL tool each replaces.
 */
#ifndef _HAL3CLI_H
#define _HAL3CLI_H

#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>

namespace hal3 {

/* A command-line mistake: reported with a pointer to --help, exit status 1. */
struct CliError : public std::runtime_error {
    explicit CliError(const std::string &what) : std::runtime_error(what) {}
};

/* Whole numbers and decimals from option values; CliError names the option. */
int64_t parseInt(const char *option, const char *value);
double parseDouble(const char *option, const char *value);

/* Run a command: exceptions become "hal <command>: message" on stderr and status 1. */
int runCommand(const char *command, const std::function<int()> &body);

/* The commands. */
int main_convert(int argc, char **argv);
int main_validate(int argc, char **argv);
int main_liftover(int argc, char **argv);
int main_maf(int argc, char **argv);
int main_fasta(int argc, char **argv);
int main_stats(int argc, char **argv);
int main_coverage(int argc, char **argv);
}

#endif
