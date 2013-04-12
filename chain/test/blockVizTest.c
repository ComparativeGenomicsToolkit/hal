
#include <stdlib.h>
#include <stdio.h>

#include "halBlockViz.h"

#ifdef ENABLE_UDC
#ifdef __cplusplus
extern "C" {
#endif
#include <pthread.h>
#include "common.h"
#include "udc.h"
#ifdef __cplusplus
}
#endif
#endif

struct bv_args_t
{
   char* path; 
   char* qSpecies; 
   char* tSpecies;
   char* tChrom; 
   int tStart; 
   int tEnd;
   int doSeq;
   char* udcPath;
};

static int parseArgs(int argc, char** argv, bv_args_t* args)
{
  if (argc != 7 && argc != 8 && argc != 9)
  {
    return -1;
  }
  args->path = argv[1];
  args->qSpecies = argv[2];
  args->tSpecies = argv[3];
  args->tChrom = argv[4];
  if (sscanf(argv[5], "%d", &args->tStart) != 1 || 
      sscanf(argv[6], "%d", &args->tEnd) != 1)
  {
    return -1;
  }
  args->doSeq = 0;
  if (argc >= 8)
  {
    if (sscanf(argv[7], "%d", &args->doSeq) != 1)
    {
      return -1;
    }
  }
  args->udcPath = NULL;
  if (argc >= 9)
  {
    args->udcPath = argv[8];
  }
  return 0; 
}

static void printBlock(FILE* file, struct hal_block_t* b)
{
  fprintf(file, "chr:%s, tSt:%d, qSt:%d, size:%d, strand:%c: %s\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand, b->sequence);
}

static void printStats(FILE* file, int handle)
{
  hal_species_t* species = halGetSpecies(handle);
  for (; species != NULL; species = species->next)
  {
    fprintf(file, "species:%s, len=%u, nc=%u, par=%s, bl=%lf\n",
            species->name, species->length, species->numChroms,
            species->parentName, species->parentBranchLength);
    hal_chromosome_t* chrom = halGetChroms(handle, species->name);
    for (; chrom != NULL; chrom = chrom->next)
    {
      int len = chrom->length > 10 ? 10 : chrom->length;
      char* dna = halGetDna(handle, species->name, chrom->name,
                            0, len);
      fprintf(file, "  chrom:%s, len=%u seq=%s...\n",
              chrom->name, chrom->length, dna);
      free(dna);              
    }
  }
}

#ifdef ENABLE_UDC
void* getBlocksWrapper(void* voidArgs)
{
  bv_args_t* args = (bv_args_t*)voidArgs;
  int handle = halOpen(args->path);
  hal_block_t* head = NULL;
  if (handle >= 0)
  {
    head = halGetBlocksInTargetRange(handle,
                                     args->qSpecies,
                                     args->tSpecies,
                                     args->tChrom, 
                                     args->tStart,
                                     args->tEnd, 
                                     args->doSeq, 
                                     0);
    halFreeBlocks(head);
  }
  pthread_exit(NULL);
}
#endif

int main(int argc, char** argv)
{
  bv_args_t args;
  
  if (parseArgs(argc, argv, &args) != 0)
  {
    fprintf(stderr, "Usage: %s <halPath> <qSpecies> <tSpecies> <tChrom> "
            "<tStart> <tEnd> [doSeq=0] [udcPath=NULL]\n\n", argv[0]);
    return -1;
  }
#ifdef ENABLE_UDC
  if (args.udcPath != NULL)
  {
    udcSetDefaultDir(args.udcPath);
  }
#endif
     
  int handle = halOpen(args.path);
  int ret = 0;
  if (handle >= 0)
  {
    printStats(stdout, handle);

    struct hal_block_t* head = halGetBlocksInTargetRange(handle, 
                                                         args.qSpecies,
                                                         args.tSpecies,
                                                         args.tChrom, 
                                                         args.tStart,
                                                         args.tEnd, 
                                                         args.doSeq, 
                                                         0);
    if (head == NULL)
    {
      ret = -1;
    }
    struct hal_block_t* cur = head;
    while (cur)
    {
      printBlock(stdout, cur);
      cur = cur->next;
    }
    halFreeBlocks(head);

#ifdef ENABLE_UDC
    #define NUM_THREADS 10
    pthread_t threads[NUM_THREADS];
    printf("\nTesting %d threads\n", NUM_THREADS);
    for (size_t t = 0; t < NUM_THREADS; ++t)
    {
      pthread_create(&threads[t], NULL, getBlocksWrapper, (void *)&args);
    }
#endif

    halClose(handle);
  }
  else
  {
    ret = -1;
  }
  return ret;
}
