
#include <stdlib.h>
#include <stdio.h>

#include "halBlockViz.h"

static int parseArgs(int argc, char** argv, char** path, char** qSpecies, 
                     char** tChrom, int* tStart, int* tEnd)
{
  if (argc != 6)
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
  return 0; 
}

static void printBlock(FILE* file, struct block* b)
{
  fprintf(file, "chr:%s, tSt:%d, qSt:%d, size:%d, strand:%c\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand);
}

int main(int argc, char** argv)
{
  char* path;
  char* qSpecies;
  char* tChrom;
  int tStart;
  int tEnd;
  
  if (parseArgs(argc, argv, &path, &qSpecies, &tChrom, &tStart, &tEnd) != 0)
  {
    fprintf(stderr, "Usage: %s <halPath> <qSpecies> <tChrom> <tStart> "
            "<tEnd>\n\n", argv[0]);
    return -1;
  }

  int handle = halOpen(path, qSpecies);
  if (handle >= 0)
  {
    struct block* head = halGetBlocksInTargetRange(handle, tChrom, tStart,
                                                   tEnd, 0);
    
    while (head)
    {
      printBlock(stdout, head);
      head = head->next;
    }

    halClose(handle);
    return 0;
  }
  
  return -1;
}
