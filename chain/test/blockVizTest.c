
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
   char* coalescenceLimit;
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
  args->coalescenceLimit = NULL;
  if (argc >= 10)
  {
    args->coalescenceLimit = argv[9];
  }
  args->udcPath = NULL;
  if (argc >= 11)
  {
    args->udcPath = argv[10];
  }
  return 0;
}

static void printBlock(FILE* file, struct hal_block_t* b)
{
  fprintf(file, "chr:%s, tSt:%ld, qSt:%ld, size:%ld, strand:%c: tgt : %s query: %s\n", 
          b->qChrom, b->tStart, b->qStart, b->size, b->strand, b->tSequence, b->qSequence);
}

static void printDupeList(FILE* file, struct hal_target_dupe_list_t* d)
{
  fprintf(file, "tDupe id:%ld qCrhom:%s\n", d->id, d->qChrom);
  for (hal_target_range_t* tr = d->tRange; tr; tr = tr->next)
  {
    fprintf(file, " tSt:%ld size:%ld\n", tr->tStart, tr->size);
  }
}

static void printStats(FILE* file, int handle)
{
  hal_species_t* species = halGetSpecies(handle, NULL);
  for (; species != NULL; species = species->next)
  {
    fprintf(file, "species:%s, len=%ld, nc=%ld, par=%s, bl=%lf\n",
            species->name, species->length, species->numChroms,
            species->parentName, species->parentBranchLength);
    hal_chromosome_t* chrom = halGetChroms(handle, species->name, NULL);
    for (; chrom != NULL; chrom = chrom->next)
    {
      int len = chrom->length > 10 ? 10 : chrom->length;
      char* dna = halGetDna(handle, species->name, chrom->name,
                            0, len, NULL);
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
    return halOpen(path, NULL);
  }
  int handle = halOpenLOD(path, NULL);
  hal_int_t maxQuery = halGetMaxLODQueryLength(handle, NULL);
  printf("Max LOD query: %lu\n", maxQuery);
  return handle;
}

#ifdef ENABLE_UDC
static void* getBlocksWrapper(void* voidArgs)
{
  bv_args_t* args = (bv_args_t*)voidArgs;
  int handle = openWrapper(args->path);
  hal_block_results_t* results = NULL;
  if (handle >= 0)
  {
    hal_seqmode_type_t sm = HAL_NO_SEQUENCE;
    if (args->doSeq != 0) sm = HAL_LOD0_SEQUENCE;
    results = halGetBlocksInTargetRange(handle,
                                        args->qSpecies,
                                        args->tSpecies,
                                        args->tChrom, 
                                        args->tStart,
                                        args->tEnd, 
                                        0,
                                        sm, 
                                        HAL_QUERY_AND_TARGET_DUPS,
                                        1,
                                        args->coalescenceLimit,
                                        NULL);
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
            "<tStart> <tEnd> [doSeq=0] [doDupes=0] [coalescenceLimit=NULL]"
            " [udcPath=NULL]\n\n", argv[0]);
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
    if (args.coalescenceLimit != NULL) {
      // Verify that the coalescence limit is possible.
      hal_species_t *coalescenceLimits = halGetPossibleCoalescenceLimits(handle,
                                                                         args.qSpecies,
                                                                         args.tSpecies,
                                                                         NULL);
      hal_species_t *curSpecies = coalescenceLimits;
      bool found = false;
      while (curSpecies != NULL) {
        if (strcmp(curSpecies->name, args.coalescenceLimit) == 0) {
          found = true;
          break;
        }
        curSpecies = curSpecies->next;
      }
      halFreeSpeciesList(coalescenceLimits);
      if (!found) {
        fprintf(stderr, "Coalescence limit %s is not in the set of valid "
                "coalescence limits for this mapping\n", args.coalescenceLimit);
        return -1;
      }
    }

    // printStats(stdout, handle);
    hal_seqmode_type_t sm = HAL_NO_SEQUENCE;
    if (args.doSeq != 0) sm = HAL_LOD0_SEQUENCE;
    struct hal_block_results_t* results = 
       halGetBlocksInTargetRange(handle, 
                                 args.qSpecies,
                                 args.tSpecies,
                                 args.tChrom, 
                                 args.tStart,
                                 args.tEnd, 
                                 0,
                                 sm, 
                                 HAL_QUERY_AND_TARGET_DUPS,
                                 1,
                                 args.coalescenceLimit,
                                 NULL);
    if (results == NULL)
    {
      fprintf(stderr, "halGetBlocksInTargetRange returned NULL\n");
      ret = -1;
    }
    else
    {
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
    }
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
