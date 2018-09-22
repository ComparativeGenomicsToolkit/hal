#include "psl_merger.h"


// Assumes a.start < b.start
bool are_syntenic(const PslBlock& a, const PslBlock& b){
    return a.qEnd <= b.qStart &&
            a.tEnd <= b.tStart &&
             a.tName == b.tName &&
                a.strand == b.strand;
}

bool is_not_overlapping_ordered_pair(const PslBlock& a, const PslBlock& b, 
                                                const hal_size_t threshold) {
    return are_syntenic(a, b) &&
             0 <= b.qStart - a.qEnd && 
                b.qStart - a.qEnd < threshold && 
                 0 <= b.tStart - a.tEnd &&
                    b.tStart - a.tEnd < threshold;
}

std::vector<int> get_next(const int pos, 
                const std::vector<PslBlock>& queryGroup,
                 const hal_size_t maxAnchorDistance) {
    std::vector<int> f;
    for (auto i = pos + 1; i < (int) queryGroup.size(); ++i) {
        if (is_not_overlapping_ordered_pair(queryGroup[pos], queryGroup[i], 
                maxAnchorDistance)) {
            if (f.empty()) f.push_back(i);
            else {
                if (not f.empty() && 
                        is_not_overlapping_ordered_pair(queryGroup[f[0]], 
                        queryGroup[i], maxAnchorDistance)) {
                    return f;
                }
                else {
                    f.push_back(i);
                }
            }
        }
    }
    
    return f;
}


// *dag* is represented as a map:  vertex -> its possible next vertices;
// *hidden_vertices* is a set of vertices that are already in paths;
// Weight of an edge equals length of the next psl block;
// Weight of a vertex equals estimated weight:
// w_j < w_i + w_e(ij) => w_j must be updated
// Also keeps track of how we came to this state:
// (prev_vertex, weight)
std::map<int, std::pair<int, hal_size_t> > 
                                weigh_dag(const std::vector<PslBlock>& group, 
                                        std::map<int, std::vector<int> >& dag, 
                                        const std::set<int>& hiddenVertices,
                                        const hal_size_t maxAnchorDistance){
    std::map<int,std::pair<int,hal_size_t> > weightedDag;
    for (int i = 0; i < (int) group.size(); ++i) {
        if (hiddenVertices.count(i)) continue; 
        std::vector<int> nexts;
        if (not dag.count(i)) {
            nexts = get_next(i, group, maxAnchorDistance);
            dag[i] = nexts;
        }
        else {
            nexts = dag[i];
       }
       // If this vertex was never visited then its weight equals to its size.
       // Intuition: otherwise we will never count its size
       if (not weightedDag.count(i))
            weightedDag[i] = std::make_pair(-1, group[i].size);
       for (auto j: nexts){
            if (hiddenVertices.count(j)) continue;
            auto alternativeWeight = weightedDag[i].second + group[j].size;
            if (not weightedDag.count(j) or 
                    weightedDag[j].second < alternativeWeight)
                // Remembers previous vertex itself and 
                // its corresponding weight:
                // previous vertex w_i + weight of the next edge 
                weightedDag[j] = std::make_pair(i, alternativeWeight);
       }
    }
    return weightedDag; 
}

int get_maxed_vertex(const std::map<int, std::pair<int, hal_size_t> >& weightedDag) {
    hal_size_t maxWeight = weightedDag.cbegin()->second.second;
    int maxPos = weightedDag.cbegin()->first;
    for (auto it = weightedDag.cbegin(); it != weightedDag.cend(); ++it) {
        if (it->second.second >= maxWeight) {
            maxWeight = it->second.second;
            maxPos = it->first;
        }
    }
    return maxPos;
}

std::vector<PslBlock> traceback(std::map<int, std::pair<int, hal_size_t> >& weightedDag, 
                                std::set<int>& hiddenVertices, 
                                const std::vector<PslBlock>& group) {
    // Chooses the path of the heaviest weight
    auto startVertex = get_maxed_vertex(weightedDag);
    std::vector<int> path = {startVertex};
    auto prevVertex = weightedDag[startVertex].first;
    while (prevVertex != -1) {
        path.push_back(prevVertex);
        prevVertex = weightedDag[prevVertex].first;
    }
    hiddenVertices.insert(path.begin(), path.end());
    std::vector<PslBlock> pslBlockPath;
    std::transform(path.begin(), path.end(), std::back_inserter(pslBlockPath), 
            [group] (int x) -> PslBlock {return group[x];});
    std::reverse(pslBlockPath.begin(),pslBlockPath.end());
    return pslBlockPath;
}
struct {    
    bool operator()(PslBlock a, PslBlock b) const
    {
        if (a.qStart < b.qStart )
            return true;
        else if (a.qStart == b.qStart) {
            return a.tStart <= b.tStart;
        }
        return false;
    }
} qStartLess;

std::vector<std::vector<PslBlock> > dag_merge(const std::vector<PslBlock>& blocks, 
                                        const hal_size_t minBlockBreath, 
                                        const hal_size_t maxAnchorDistance){
    std::map<std::string, std::vector<PslBlock> > blocksByQName;
    for (auto block : blocks) blocksByQName[block.qName].push_back(block);
    std::vector<std::vector<PslBlock> > paths;
    for (auto pairs : blocksByQName) {
        std::vector<PslBlock> group = pairs.second; 
        std::string qName = pairs.first;
        std::map<int, std::vector<int> > dag;
        std::set<int> hiddenVertices;
        std::sort(group.begin(), group.end(), qStartLess); 
        while (hiddenVertices.size() != group.size()){
            auto weightedDag = weigh_dag(group, dag, hiddenVertices, 
                    maxAnchorDistance);
            auto path = traceback(weightedDag, hiddenVertices, group);
            if (path.empty()) {
                break;
            }
            auto qLen = path.back().qEnd - path[0].qStart;
            auto tLen = path.back().tEnd - path[0].tStart;
            if (qLen >= minBlockBreath && tLen >= minBlockBreath){
            //if (qLen >= minBlockBreath && tLen >= minBlockBreath){
                paths.push_back(path);
            }
        }
    }
    return paths;
   
}



