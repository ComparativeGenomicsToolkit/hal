
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
   int doDupes;
   char* udcPath;
};

static int parseArgs(int argc, char** argv, bv_args_t* args)
{
  if (argc != 7 && argc != 8 && argc != 9 && argc != 10)
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
  args->doDupes = 0;
  if (argc >= 9)
  {
    if (sscanf(argv[8], "%d", &args->doDupes) != 1)
    {
      return -1;
    }
  }
  args->udcPath = NULL;
  if (argc >= 10)
  {
    args->udcPath = argv[9];
  }
  return 0; 
}

static void printBlock(FILE* file, struct hal_block_t* b)
{
  fprintf(file, "chr:%s, tSt:%ld, qSt:%ld, size:%ld, strand:%c: %s\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand, b->sequence);
}

static void printDupeList(FILE* file, struct hal_target_dupe_list_t* d)
{
  fprintf(file, "tDupe id:%d qCrhom:%s\n", d->id, d->qChrom);
  for (hal_target_range_t* tr = d->tRange; tr; tr = tr->next)
  {
    fprintf(file, " tSt:%ld size:%ld\n", tr->tStart, tr->size);
  }
}

static void printStats(FILE* file, int handle)
{
  hal_species_t* species = halGetSpecies(handle);
  for (; species != NULL; species = species->next)
  {
    fprintf(file, "species:%s, len=%ld, nc=%ld, par=%s, bl=%lf\n",
            species->name, species->length, species->numChroms,
            species->parentName, species->parentBranchLength);
    hal_chromosome_t* chrom = halGetChroms(handle, species->name);
    for (; chrom != NULL; chrom = chrom->next)
    {
      int len = chrom->length > 10 ? 10 : chrom->length;
      char* dna = halGetDna(handle, species->name, chrom->name,
                            0, len);
      fprintf(file, "  chrom:%s, len=%ld seq=%s...\n",
              chrom->name, chrom->length, dna);
      free(dna);              
    }
  }
}

static int openWrapper(char* path)
{
  if (strcmp(path + strlen(path) - 3, "hal") == 0)
  {
    return halOpen(path);
  }
  return halOpenLOD(path);
}

#ifdef ENABLE_UDC
static void* getBlocksWrapper(void* voidArgs)
{
  bv_args_t* args = (bv_args_t*)voidArgs;
  int handle = openWrapper(args->path);
  hal_block_results_t* results = NULL;
  if (handle >= 0)
  {
    results = halGetBlocksInTargetRange(handle,
                                     args->qSpecies,
                                     args->tSpecies,
                                     args->tChrom, 
                                     args->tStart,
                                     args->tEnd, 
                                     args->doSeq, 
                                     0);
    halFreeBlockResults(results);
  }
  pthread_exit(NULL);
}
#endif

int main(int argc, char** argv)
{
  bv_args_t args;
  
  if (parseArgs(argc, argv, &args) != 0)
  {
    fprintf(stderr, "Usage: %s <halLodPath> <qSpecies> <tSpecies> <tChrom> "
            "<tStart> <tEnd> [doSeq=0] [doDupes=0] [udcPath=NULL]\n\n", argv[0]);
    return -1;
  }
#ifdef ENABLE_UDC
  if (args.udcPath != NULL)
  {
    udcSetDefaultDir(args.udcPath);
  }
#endif
     
  int handle = openWrapper(args.path);
  int ret = 0;
  if (handle >= 0)
  {
    //printStats(stdout, handle);

    struct hal_block_results_t* results = 
       halGetBlocksInTargetRange(handle, 
                                 args.qSpecies,
                                 args.tSpecies,
                                 args.tChrom, 
                                 args.tStart,
                                 args.tEnd, 
                                 args.doSeq, 
                                 args.doDupes);
    if (results == NULL)
    {
      ret = -1;
    }
    struct hal_block_t* cur = results->mappedBlocks;
    while (cur)
    {
      printBlock(stdout, cur);
      cur = cur->next;
    }
    struct hal_target_dupe_list_t* dupeList = results->targetDupeBlocks;
    while (dupeList)
    {
      printDupeList(stdout, dupeList);
      dupeList = dupeList->next;
    }
    halFreeBlockResults(results);

#ifdef ENABLE_UDC
    #define NUM_THREADS 10
    pthread_t threads[NUM_THREADS];
    printf("\nTesting %d threads\n", NUM_THREADS);
    for (size_t t = 0; t < NUM_THREADS; ++t)
    {
      pthread_create(&threads[t], NULL, getBlocksWrapper, (void *)&args);
    }
#endif

    
  }
  else
  {
    ret = -1;
  }
#ifdef ENABLE_UDC
  pthread_exit(NULL);
#endif
  return ret;
}
