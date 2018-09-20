#ifndef PSL_IO_H
#define PSL_IO_H

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "psl.h"


namespace psl_io {
std::vector<std::string> split(const std::string &s, char delim);


std::vector<PslBlock> get_blocks_set(const std::string psl);

std::vector<int> get_qInserts(const std::vector<PslBlock>& blocks);

std::vector<int> get_tInserts(const std::vector<PslBlock>& blocks);

 Psl construct_psl(std::vector<PslBlock> blocks);
 
 void write_psl(const std::vector<std::vector<PslBlock> >& merged_blocks, 
         const std::string& outFilePath);
 
 }

#endif /* PSL_IO_H */


