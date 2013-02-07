
#include <stdlib.h>
#include <stdio.h>

#include "halBlockViz.h"

static int parseArgs(int argc, char** argv, char** path, char** qSpecies, 
                     char** tChrom, int* tStart, int* tEnd, int* doSeq)
{
  if (argc != 6 && argc != 7)
  {
    return -1;
  }
  *path = argv[1];
  *qSpecies = argv[2];
  *tChrom = argv[3];
  if (sscanf(argv[4], "%d", tStart) != 1 || 
      sscanf(argv[5], "%d", tEnd) != 1)
  {
    return -1;
  }
  *doSeq = 0;
  if (argc == 7)
  {
    if (sscanf(argv[6], "%d", doSeq) != 1)
    {
      return -1;
    }
  }
  return 0; 
}

static void printBlock(FILE* file, struct block* b)
{
  fprintf(file, "chr:%s, tSt:%d, qSt:%d, size:%d, strand:%c: %s\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand, b->sequence);
}

int main(int argc, char** argv)
{
  char* path;
  char* qSpecies;
  char* tChrom;
  int tStart;
  int tEnd;
  int doSeq;
  
  if (parseArgs(argc, argv, &path, &qSpecies, &tChrom, &tStart, &tEnd,
        &doSeq) != 0)
  {
    fprintf(stderr, "Usage: %s <halPath> <qSpecies> <tChrom> <tStart> "
            "<tEnd> [doSeq=0]\n\n", argv[0]);
    return -1;
  }

  int handle = halOpen(path, qSpecies);
  if (handle >= 0)
  {
    struct block* head = halGetBlocksInTargetRange(handle, tChrom, tStart,
                                                   tEnd, doSeq);
    struct block* cur = head;
    while (cur)
    {
      printBlock(stdout, cur);
      cur = cur->next;
    }
    halFreeBlocks(head);
    halClose(handle);
    return 0;
  }
  
  return -1;
}
