
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
  args->doDupes = 0;
  if (argc >= 8)
  {
    if (sscanf(argv[7], "%d", &args->doDupes) != 1)
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

static int openWrapper(char* path)
{
  if (strcmp(path + strlen(path) - 3, "hal") == 0)
  {
    return halOpen(path, NULL);
  }
  return halOpenLOD(path, NULL);
}

int main(int argc, char** argv)
{
  bv_args_t args;
  
  if (parseArgs(argc, argv, &args) != 0)
  {
    fprintf(stderr, "Usage: %s <halLodPath> <qSpecies> <tSpecies> <tChrom> "
            "<tStart> <tEnd> [doDupes=0] [udcPath=NULL]\n\n", argv[0]);
    return -1;
  }
#ifdef ENABLE_UDC
  if (args.udcPath != NULL)
  {
    udcSetDefaultDir(args.udcPath);
  }
#endif
     
  int handle = openWrapper(args.path);
  int ret = -1;
  if (handle >= 0)
  {
    // printStats(stdout, handle);
    hal_species_t qSpecies;
    qSpecies.name = args.qSpecies;
    qSpecies.next = NULL;

    long numBytes = halGetMAF(stdout,  
                              handle,
                              &qSpecies,
                              args.tSpecies,
                              args.tChrom, 
                              args.tStart,
                              args.tEnd,
                              0, // maxRefGap
                              0, // maxBlockLength
                              args.doDupes,
                              NULL);

    if (numBytes >= 0)
    {
      ret = 0;
    }
    printf("\nread %ld bytes\n", numBytes);
  }
  return ret;
}
