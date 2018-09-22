/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   hal2psl.h
 * Author: admin
 *
 * Created on July 11, 2018, 12:15 PM
 */

#ifndef HAL2PSL_H
#define HAL2PSL_H

#include "halBedLine.h"
#include "halBlockLiftover.h"

namespace hal {
class Hal2Psl : public hal::BlockLiftover {
    
    void storePslResults(std::vector<PslBlock>& pslBlocks);
    void makeUpPsl(const std::vector<hal::PSLInfo>& vpsl,
                                const std::vector<hal::BedBlock>& blocks, 
                                const char strand,
                                const hal_index_t start,
                                const std::string chrName,
                                std::vector<PslBlock>& pslBlocks);
    
    public:

    Hal2Psl(){}
    std::vector<PslBlock> convert2psl(hal::AlignmentConstPtr alignment,
                       const hal::Genome* srcGenome,
                       const hal::Genome* tgtGenome,
                       const std::string srcChrom);
};
}
#endif /* HAL_MERGER_H */

