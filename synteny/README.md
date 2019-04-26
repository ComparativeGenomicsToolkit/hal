DAG-based Reconstruction of Synteny Blocks
=====

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
* Create synteny blocks in psl format for the alignment cactus.hal including genomes Genom1 and Genome2

    `hal2synteny cactus.hal out.psl --queryGenome Genome2 --targetGenome Genome1 --maxAnchorDistance 1000000 --minBlockSize 1000000`

* Create synteny blocks in psl format just for query chromosome chrA 

    `hal2synteny cactus.hal out.psl --queryGenome Genome2 --targetGenome Genome1 --queryChromosome chrA --maxAnchorDistance 1000000 --minBlockSize 1000000`

Code Contributors
-----
* Ksenia Krasheninnikova (Theodosius Dobzhansky Center for Genome Bioinformatics)
* Joel Armstrong (UCSC)
* Mark Diekhans (UCSC)


