
#include <stdlib.h>
#include <stdio.h>

#include "halBlockViz.h"

#ifdef ENABLE_UDC
#ifdef __cplusplus
extern "C" {
#endif
#include "common.h"
#include "udc.h"
#ifdef __cplusplus
}
#endif
#endif

static int parseArgs(int argc, char** argv, char** path, char** qSpecies, 
                     char** tSpecies,
                     char** tChrom, int* tStart, int* tEnd, int* doSeq,
                     char** udcPath)
{
  if (argc != 7 && argc != 8 && argc != 9)
  {
    return -1;
  }
  *path = argv[1];
  *qSpecies = argv[2];
  *tSpecies = argv[3];
  *tChrom = argv[4];
  if (sscanf(argv[5], "%d", tStart) != 1 || 
      sscanf(argv[6], "%d", tEnd) != 1)
  {
    return -1;
  }
  *doSeq = 0;
  if (argc >= 8)
  {
    if (sscanf(argv[7], "%d", doSeq) != 1)
    {
      return -1;
    }
  }
  *udcPath = NULL;
  if (argc >= 9)
  {
    *udcPath = argv[8];
  }
  return 0; 
}

static void printBlock(FILE* file, struct hal_block_t* b)
{
  fprintf(file, "chr:%s, tSt:%d, qSt:%d, size:%d, strand:%c: %s\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand, b->sequence);
}

static void printStats(FILE* file, int handle, int threadID)
{
  hal_species_t* species = halGetSpecies(handle, threadID);
  for (; species != NULL; species = species->next)
  {
    fprintf(file, "species:%s, len=%u, nc=%u, par=%s, bl=%lf\n",
            species->name, species->length, species->numChroms,
            species->parentName, species->parentBranchLength);
    hal_chromosome_t* chrom = halGetChroms(handle, threadID, species->name);
    for (; chrom != NULL; chrom = chrom->next)
    {
      int len = chrom->length > 10 ? 10 : chrom->length;
      char* dna = halGetDna(handle, threadID, species->name, chrom->name,
                            0, len);
      fprintf(file, "  chrom:%s, len=%u seq=%s...\n",
              chrom->name, chrom->length, dna);
      free(dna);              
    }
  }
}

int main(int argc, char** argv)
{
  char* path;
  char* qSpecies;
  char* tSpecies;
  char* tChrom;
  int tStart;
  int tEnd;
  int doSeq;
  char* udcPath;
  
  if (parseArgs(argc, argv, &path, &qSpecies, &tSpecies, 
                &tChrom, &tStart, &tEnd,
                &doSeq, &udcPath) != 0)
  {
    fprintf(stderr, "Usage: %s <halPath> <qSpecies> <tSpecies> <tChrom> "
            "<tStart> <tEnd> [doSeq=0] [udcPath=NULL]\n\n", argv[0]);
    return -1;
  }
#ifdef ENABLE_UDC
  if (udcPath != NULL)
  {
    udcSetDefaultDir(udcPath);
  }
#endif
  int ret = halInit(1);
  int handle = halOpen(path, 0);
  if (ret == 0 && handle >= 0)
  {
    printStats(stdout, handle, 0);

    struct hal_block_t* head = halGetBlocksInTargetRange(handle, 
                                                         0,
                                                         qSpecies,
                                                         tSpecies,
                                                         tChrom, 
                                                         tStart,
                                                         tEnd, 
                                                         doSeq, 
                                                         0);
    struct hal_block_t* cur = head;
    while (cur)
    {
      printBlock(stdout, cur);
      cur = cur->next;
    }
    halFreeBlocks(head);
    halClose(handle, 0);
  }
  else
  {
    ret = -1;
  }
  halExit();
  return ret;
}
