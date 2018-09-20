/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   psl_merger.h
 * Author: admin
 *
 * Created on July 6, 2018, 1:13 PM
 */

#ifndef PSL_MERGER_H
#define PSL_MERGER_H

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <ctime> 

#include "psl.h"
bool are_syntenic(const PslBlock& a, const PslBlock& b);

bool is_not_overlapping_ordered_pair(const PslBlock& a, const PslBlock& b, 
                                                const hal_size_t threshold=5000);

std::vector<int> get_next(const int pos, 
                const std::vector<PslBlock>& queryGroup,
                 const hal_size_t maxAnchorDistance=5000);


std::map<int, std::pair<int, hal_size_t> > 
                                weigh_dag(const std::vector<PslBlock>& group, 
                                        std::map<int, std::vector<int> >& dag, 
                                        const std::set<int>& hiddenVertices,
                                        const hal_size_t maxAnchorDistance);

int get_maxed_vertex(const std::map<int, std::pair<int, int> >& weightedDag);

std::vector<PslBlock> traceback(std::map<int, std::pair<int, hal_size_t> >& weightedDag, 
                                std::set<int>& hiddenVertices, 
                                const std::vector<PslBlock>& group) ;


std::vector<std::vector<PslBlock> > 
                                dag_merge(const std::vector<PslBlock>& blocks, 
                                        const hal_size_t minBlockBreath, 
                                        const hal_size_t maxAnchorDistance);



#endif /* PSL_MERGER_H */

