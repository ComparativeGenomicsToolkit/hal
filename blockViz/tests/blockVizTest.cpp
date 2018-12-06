/*
 * This is done mostly using C style since testing C interface, done in C++ to
 * use option parser.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "halBlockViz.h"
#include "halCLParser.h"
#include "halAlignmentInstance.h"

// for debugging
#undef UDC_DEBUG_VERBOSE

#ifdef ENABLE_UDC
#ifdef __cplusplus
extern "C" {
#endif
#include <pthread.h>
#include "common.h"
#include "udc.h"
#ifdef UDC_DEBUG_VERBOSE
#include "verbose.h"
#endif
#ifdef __cplusplus
}
#endif
#endif

struct bv_args_t {
   char* path; 
   char* qSpecies; 
   char* tSpecies;
   char* tChrom; 
   int tStart; 
   int tEnd;
   int doSeq;
   int doDupes;
   char* udcCacheDir;
   char* coalescenceLimit;
    int verbose;
};

static void initParser(hal::CLParser& optionsParser) {
    optionsParser.setDescription("Test blockViz code from command line");
    optionsParser.addOptionFlag("verbose", "verbose tracing", false);
    optionsParser.addOptionFlag("doSeq", "get seqeuence", false);
    optionsParser.addOptionFlag("doDupes", "get duplicate regions",false);
    optionsParser.addOption("coalescenceLimit", " coalescence limit specices, default is none", "");
    optionsParser.addArgument("halLodPath", "path to HAL or LOD file");
    optionsParser.addArgument("qSpecies", "query species name");
    optionsParser.addArgument("tSpecies", "target species name");
    optionsParser.addArgument("tChrom", "target chromosome");
    optionsParser.addArgument("tStart", "zero based start in target");
    optionsParser.addArgument("tEnd", "half-open end in target");
}

static char* argStr(hal::CLParser& optionsParser,
                    const std::string& name) {
    // need to leak memory, as string is transient
    std::string* val = new std::string(optionsParser.getArgument<std::string>(name));
    return const_cast<char*>(val->c_str());
}

static char* optionStrOrNull(hal::CLParser& optionsParser,
                             const std::string& name) {
    // need to leak memory, as string is transient
    std::string* val = new std::string(optionsParser.getOption<std::string>(name));
    if (val->empty()) {
        return NULL;
    } else {
        return const_cast<char*>(val->c_str());
    }
}

static bool parseArgs(int argc, char** argv, bv_args_t* args) {
    hal::CLParser optionsParser(hal::READ_ACCESS);
    initParser(optionsParser);
    try {
        optionsParser.parseOptions(argc, argv);
    } catch (hal_exception &e) {
        std::cerr << e.what() << std::endl;
        optionsParser.printUsage(std::cerr);
        return false;
    }
    args->path = argStr(optionsParser, "halLodPath");
    args->qSpecies = argStr(optionsParser, "qSpecies");
    args->tSpecies = argStr(optionsParser, "tSpecies");
    args->tChrom = argStr(optionsParser, "tChrom");
    args->tStart = optionsParser.get<int>("tStart");
    args->tEnd = optionsParser.get<int>("tEnd");
    args->doSeq = optionsParser.get<bool>("doSeq");
    args->doDupes = optionsParser.get<bool>("doDupes") ;
    args->udcCacheDir = optionStrOrNull(optionsParser, "udcCacheDir");
    args->coalescenceLimit = optionStrOrNull(optionsParser, "coalescenceLimit");
    args->verbose = optionsParser.get<bool>("verbose") ;
    return true;
}

static void printBlock(FILE* file, struct hal_block_t* b)
{
  fprintf(file, "chr:%s, tSt:%ld, qSt:%ld, size:%ld, strand:%c: tgt : %.10s query: %.10s\n", 
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

static int openWrapper(char* path)
{
    if (not hal::detectHalAlignmentFormat(path).empty()) {
        return halOpen(path, NULL);
    } else {
        int handle = halOpenLOD(path, NULL);
        hal_int_t maxQuery = halGetMaxLODQueryLength(handle, NULL);
        printf("Max LOD query: %lu\n", maxQuery);
        return handle;
    }
}

/* Verify that the coalescence limit is possible. */
static bool checkCoalescenceLimit(int handle,
                                  struct bv_args_t* args) {
      hal_species_t *coalescenceLimits = halGetPossibleCoalescenceLimits(handle,
                                                                         args->qSpecies,
                                                                         args->tSpecies,
                                                                         NULL);
      hal_species_t *curSpecies = coalescenceLimits;
      bool found = false;
      while ((curSpecies != NULL) && !found)  {
        if (strcmp(curSpecies->name, args->coalescenceLimit) == 0) {
          found = true;
        } else {
            curSpecies = curSpecies->next;
        }
      }
      halFreeSpeciesList(coalescenceLimits);
      if (!found) {
          fprintf(stderr, "Coalescence limit %s is not in the set of valid "
                  "coalescence limits for this mapping\n", args->coalescenceLimit);
      }
      return found;
}

#ifdef ENABLE_UDC
static bool someThreadFailed = false;

static void* getBlocksWrapper(void* voidArgs)
{
    struct bv_args_t* args = static_cast<struct bv_args_t*>(voidArgs);
  int handle = openWrapper(args->path);
  if (handle < 0) {
      someThreadFailed = true;
  }  else {
      hal_block_results_t* results = NULL;
      hal_seqmode_type_t sm = HAL_NO_SEQUENCE;
      if (args->doSeq != 0)
          sm = HAL_LOD0_SEQUENCE;
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
      if (results == NULL) {
          someThreadFailed = true;
      }
      halFreeBlockResults(results);
  }
  pthread_exit(NULL);
}
#endif

static bool runTest(bv_args_t* args,
                    int handle) {
    if (args->coalescenceLimit != NULL) {
        if (!checkCoalescenceLimit(handle, args)) {
            return false;
        }
    }
    hal_seqmode_type_t sm = HAL_NO_SEQUENCE;
    if (args->doSeq != 0) {
        sm = HAL_LOD0_SEQUENCE;
    }
    struct hal_block_results_t* results = 
        halGetBlocksInTargetRange(handle, 
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
    if (results == NULL)
    {
        fprintf(stderr, "halGetBlocksInTargetRange returned NULL\n");
        return false;
    }
    if (args->verbose) {
    struct hal_block_t* cur = results->mappedBlocks;
        while (cur != NULL)
            {
                printBlock(stdout, cur);
                cur = cur->next;
            }
        struct hal_target_dupe_list_t* dupeList = results->targetDupeBlocks;
        while (dupeList != NULL)
            {
                printDupeList(stdout, dupeList);
                dupeList = dupeList->next;
            }
    }
    halFreeBlockResults(results);
#ifdef ENABLE_UDC
    #define NUM_THREADS 10
    pthread_t threads[NUM_THREADS];
    printf("\nTesting %d threads\n", NUM_THREADS);
    for (size_t t = 0; t < NUM_THREADS; ++t)
    {
      pthread_create(&threads[t], NULL, getBlocksWrapper, args);
    }
  pthread_exit(NULL);
  if (someThreadFailed) {
      return false;
  }
#endif
    return true;
}

int main(int argc, char** argv)
{
  bv_args_t args;
  if (!parseArgs(argc, argv, &args)) {
      return 1;
  }
    
#ifdef ENABLE_UDC
#ifdef UDC_DEBUG_VERBOSE
  verboseSetLevel(100);
#endif
    if (args.udcCacheDir)
  {
      udcSetDefaultDir(args.udcCacheDir);
  }
#endif
    int handle = openWrapper(args.path);
    if (handle < 0) {
        return 1;
    }
    if (!runTest(&args, handle)) {
        return 1;
    }
  return 0;
}
