DAG-based Reconstruction of Synteny Blocks
=====

Run-time parameters
-----
 Option |  Effect
--- | ---
`--maxAnchorDistance <value>`  | upper bound on distance for syntenic blocks, default is 5Kb 
`--minBlockSize <value>`        | lower bound on synteny block length, default is 5Kb 
`--queryChromosome <value>`     | chromosome to infer synteny, default is whole genome 
`--queryGenome <value>`         | source genome name 
`--targetGenome <value>`        | reference genome name 

Other parameters describe the options for handling alignment file, in particular HDF5 options (applicable if alignments are in HAL format). Detailed information can be found in the main [README.md](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md)
    
Algorithm
-----
In order to build a set of most continuous synteny blocks covering as much of genomes as possible we build a graph and apply the algorithm:
    
1. Initialize weight labels of vertices and edges:
    
* initial weight of each vertex is initialized as the size of the corresponding alignment block
    
* the weight of each edge coming into a vertex i equals to the length of the corresponding alignment block (initial weight of the vertex i)
    
* for each vertex A consider its possible candidate blocks (descendants), update weight of each descendant vertex B in case the weight of an edge coming into B + weight of the previous vertex A is greater than weight of B; if the weight was changed update the id of the previous vertex
    
2. Find the vertex with maximal weight and trace back
    
3. Remove the vertices that comprise the best path from the graph
    
4. If not all vertices are in some paths then go to 2

Sample Usage
-----
* Create synteny blocks for the alignment cactus.hal including genomes Genome1 and Genome2

    `halSynteny  --queryGenome Genome2 --targetGenome Genome1 --maxAnchorDistance 1000000 --minBlockSize 1000000 cactus.hal out.psl`

* Create synteny blocks for alignments in [PSL](http://genome.ucsc.edu/FAQ/FAQformat#format2) format from for query chromosome chrA 

    `halSynteny --queryGenome Genome2 --targetGenome Genome1 --queryChromosome chrA --maxAnchorDistance 1000000 --minBlockSize 1000000 cactus.hal out.psl`

Code Contributors
-----
* Ksenia Krasheninnikova
* Joel Armstrong 
* Mark Diekhans 


