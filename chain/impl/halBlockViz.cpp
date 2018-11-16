/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <map>
#include <sstream>
#include <cmath>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <limits>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hal.h"
#include "halChain.h"
#include "halBlockViz.h"
#include "halBlockMapper.h"
#include "halLodManager.h"
#include "halMafExport.h"

#ifdef ENABLE_UDC
#include <pthread.h>
static pthread_mutex_t HAL_MUTEX = PTHREAD_MUTEX_INITIALIZER;
static inline void halLock() {
    if (pthread_mutex_lock(&HAL_MUTEX) != 0) {
        throw hal_exception("pthread_mutex_lock failed: " + std::string(std::strerror(errno)));
    }
}
static inline void halUnlock() {
    if (pthread_mutex_unlock(&HAL_MUTEX) != 0) {
        throw hal_exception("pthread_mutex_unlock failed: " + std::string(std::strerror(errno)));
    }
}
#else
static inline void halLock() {
}
static inline void halUnlock() {
}
#endif

using namespace std;
using namespace hal;

typedef map<int, pair<string, LodManagerPtr> > HandleMap;
static HandleMap handleMap;

static int halOpenLodOrHal(char* inputPath, bool isLod, char **errStr);
static void checkHandle(int handle);
static void checkGenomes(int halHandle, 
                         const Alignment* alignment, const string& qSpecies,
                         const string& tSpecies, const string& tChrom);

static const Alignment* getExistingAlignment(int handle,
                                              hal_size_t queryLength,
                                              bool needSequence);
static bool isAlignmentLod0(int handle, hal_size_t queryLength);
static char* copyCString(const string& inString);

static hal_block_results_t* readBlocks(const Alignment* seqAlignment,
                                       const Sequence* tSequence,
                                       hal_index_t absStart, 
                                       hal_index_t absEnd,
                                       bool tReversed,
                                       const Genome* qGenome, 
                                       bool getSequenceString,
                                       bool doDupes, bool doTargetDupes,
                                       bool doAdjes, const char *coalescenceLimitName);

static void readBlock(const Alignment* seqAlignment,
                      hal_block_t* cur, 
                      vector<MappedSegmentPtr>& fragments,                                        bool getSequenceString, const string& genomeName);

static hal_target_dupe_list_t* processTargetDupes(BlockMapper& blockMapper,
                                                  MappedSegmentSet& paraSet);

static void readTargetRange(hal_target_dupe_list_t* cur,
                            vector<MappedSegmentPtr>& fragments);

static void cleanTargetDupesList(vector<hal_target_dupe_list_t*>& dupeList);

static void mergeCompatibleDupes(vector<hal_target_dupe_list_t*>& dupeList);


/* return an error in errStr if not NULL, otherwise throw and exception with
 * the message. */
static void handleError(const string& msg,
                        char **errStr) {
    if (errStr == NULL) {
        throw hal_exception(msg);
    } else {
        *errStr = stString_copy(msg.c_str());
    }
}

extern "C" int halOpenLOD(char *lodFilePath, char **errStr)
{
  bool isHal = lodFilePath && strlen(lodFilePath) > 4 &&
     strcmp(lodFilePath + strlen(lodFilePath) - 4, ".hal") == 0;

  return halOpenLodOrHal(lodFilePath, !isHal, errStr);
}

extern "C" int halOpen(char* halFilePath, char **errStr)
{
  return halOpenLodOrHal(halFilePath, false, errStr);
}

int halOpenLodOrHal(char* inputPath, bool isLod, char **errStr)
{
  halLock();
  int handle = -1;
  try
  {
    for (HandleMap::iterator mapIt = handleMap.begin(); 
         mapIt != handleMap.end(); ++mapIt)
    {
      if (mapIt->second.first == string(inputPath))
      {
        handle = mapIt->first;
      }
    }
    if (handle == -1)
    {
      HandleMap::reverse_iterator mapIt = handleMap.rbegin();
      if (mapIt == handleMap.rend())
      {
        handle = 0;
      }
      else
      {
        handle = mapIt->first + 1;
      }
      LodManagerPtr lodManager(new LodManager());
      if (isLod == true)
      {
        lodManager->loadLODFile(inputPath);
      }
      else
      {
        lodManager->loadSingeHALFile(inputPath);
      }
      handleMap.insert(pair<int, pair<string, LodManagerPtr> >(
                         handle, pair<string, LodManagerPtr>(
                           inputPath, lodManager)));
    }
  }
  catch(exception& e)
  {
      halUnlock();
      handleError("halOpenLodOrHal error: " + string(inputPath) + ": " + e.what(), errStr);
      return -1;
  }
  catch(...)
  {
      halUnlock();
    handleError("halOpenLodOrHal error: " + string(inputPath) + ": Unknown exception", errStr);
    return -1;
  }
  halUnlock();
  return handle;
}

extern "C" int halClose(int handle, char **errStr)
{
  halLock();
  int ret = 0;
  try
  {
    HandleMap::iterator mapIt = handleMap.find(handle);
    if (mapIt == handleMap.end())
    {
        halUnlock();
        handleError("halClose error on handle: " + std::to_string(handle) + ": not found", errStr);
        return -1;
    }
    handleMap.erase(mapIt);
  }
  catch(exception& e)
  {
      halUnlock();
      handleError("halClose error on handle: " + std::to_string(handle) + ": " + e.what(), errStr);
      return -1;
  }
  catch(...)
  {
    halUnlock();
    handleError("halClose error on handle: " + std::to_string(handle) + ": unknown exception", errStr);
    return -1;
  }
  halUnlock();
  return ret;
}

extern "C" void halFreeBlockResults(struct hal_block_results_t* results)
{
  if (results != NULL)
  {
    halFreeBlocks(results->mappedBlocks);
    halFreeTargetDupeLists(results->targetDupeBlocks);
    free(results);
  }
}

extern "C" void halFreeBlocks(struct hal_block_t* head)
{
  while (head != NULL)
  {
    hal_block_t* next = head->next;
    free(head->qSequence);
    free(head->tSequence);
    free(head->qChrom);
    free(head);
    head = next;
  }
}

extern "C" void halFreeTargetDupeLists(struct hal_target_dupe_list_t* dupes)
{
  while (dupes != NULL)
  {
    hal_target_dupe_list_t* next = dupes->next;
    while (dupes->tRange != NULL)
    {
      hal_target_range_t* nextRange = dupes->tRange->next;
      free(dupes->tRange);
      dupes->tRange = nextRange;
    }
    free(dupes->qChrom);
    free(dupes);
    dupes = next;
  }
}

extern "C" 
struct hal_block_results_t *halGetBlocksInTargetRange(int halHandle,
                                                      char* qSpecies,
                                                      char* tSpecies,
                                                      char* tChrom,
                                                      hal_int_t tStart, 
                                                      hal_int_t tEnd,
                                                      hal_int_t tReversed,
                                                      hal_seqmode_type_t seqMode,
                                                      hal_dup_type_t dupMode,
                                                      int mapBackAdjacencies,
                                                      const char *coalescenceLimitName,
                                                      char **errStr)
{
  halLock();
  hal_block_results_t* results = NULL;
  try
  {
    hal_int_t rangeLength = tEnd - tStart;
    if (rangeLength < 0)
    {
        halUnlock();
        handleError("halGetBlocksInTargetRange invalid query range [" + std::to_string(tStart) + "," + std::to_string(tEnd) + ")", errStr);
        return NULL;
    }
    if (tReversed != 0 && mapBackAdjacencies != 0)
    {
      halUnlock();
      handleError("halGetBlocksInTargetRange tReversed can only be set when mapBackAdjacencies is 0", errStr);
      return NULL;
    }
    if (tReversed != 0 && dupMode == HAL_QUERY_AND_TARGET_DUPS)
    {
        halUnlock();
        handleError("tReversed cannot be set in conjunction with dupMode=HAL_QUERY_AND_TARGET_DUPS", errStr);
        return NULL;
    }
    bool getSequenceString;
    switch (seqMode) 
    {
    case HAL_NO_SEQUENCE: getSequenceString = false; break;
    case HAL_FORCE_LOD0_SEQUENCE: getSequenceString = true; break;           
    case HAL_LOD0_SEQUENCE: default:
      getSequenceString = isAlignmentLod0(halHandle, hal_size_t(rangeLength)); 
    }
      
    const Alignment* alignment = getExistingAlignment(halHandle, hal_size_t(rangeLength), getSequenceString);
    checkGenomes(halHandle, alignment, qSpecies, tSpecies, tChrom);

    const Genome* qGenome = alignment->openGenome(qSpecies);
    const Genome* tGenome = alignment->openGenome(tSpecies);
    const Sequence* tSequence = tGenome->getSequence(tChrom);

    hal_index_t myEnd = tEnd > 0 ? tEnd : tSequence->getSequenceLength();
    hal_index_t absStart = tSequence->getStartPosition() + tStart;
    hal_index_t absEnd = tSequence->getStartPosition() + myEnd - 1;
    if (absStart > absEnd)
    {
        halUnlock();
        handleError("halGetBlocksInTargetRange invalid range", errStr);
        return NULL;
    }
    if (absEnd > tSequence->getEndPosition())
    {
        halUnlock();
        handleError("halGetBlocksInTargetRange target end position outside of target sequence", errStr);
        return NULL;
    }
    // We now know the query length so we can do a proper lod query
    if (tEnd == 0)
    {
        alignment = getExistingAlignment(halHandle, absEnd - absStart, false);
        checkGenomes(halHandle, alignment, qSpecies, tSpecies, tChrom);
      qGenome = alignment->openGenome(qSpecies);
      tGenome = alignment->openGenome(tSpecies);
      tSequence = tGenome->getSequence(tSequence->getName());
    }

    const Alignment* seqAlignment = NULL;
    if (getSequenceString == true)
    {
      // note: this separate pointer no longer necessary since we will
      // not get sequence unless alignment has sequence.  don't bother
      // getting rid of it since it allows us to easily revert back to 
      // the previous functionaly of allowing lod-blocks to acces lod-0
      // sequence
        seqAlignment = getExistingAlignment(halHandle, absEnd - absStart, true);
    }

    results = readBlocks(seqAlignment, tSequence, absStart, absEnd, 
                         tReversed != 0,
                         qGenome,
                         getSequenceString, dupMode != HAL_NO_DUPS, 
                         dupMode == HAL_QUERY_AND_TARGET_DUPS,
                         mapBackAdjacencies != 0, coalescenceLimitName);
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetBlocksInTargetRange error reading blocks: " + string(e.what()), errStr);
    return NULL;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetBlocksInTargetRange error reading blocks: unknown exception" , errStr);
    return NULL;
  }
  halUnlock();
  return results;
}

extern "C"
struct hal_block_results_t *halGetBlocksInTargetRange_filterByChrom(int halHandle,
                                                                    char* qSpecies,
                                                                    char* tSpecies,
                                                                    char* tChrom,
                                                                    hal_int_t tStart,
                                                                    hal_int_t tEnd,
                                                                    hal_int_t tReversed,
                                                                    hal_seqmode_type_t seqMode,
                                                                    hal_dup_type_t dupMode,
                                                                    int mapBackAdjacencies,
                                                                    char *qChrom,
                                                                    const char *coalescenceLimitName,
                                                                    char **errStr) {
    hal_block_results_t* results = halGetBlocksInTargetRange(halHandle, qSpecies,
                                                             tSpecies, tChrom,
                                                             tStart, tEnd,
                                                             tReversed, seqMode,
                                                             dupMode, mapBackAdjacencies,
                                                             coalescenceLimitName, errStr);
    if (results == NULL) {
        // An error occured.
        return NULL;
    }

    // Filter out any blocks not on qChrom.
    hal_block_t *blockHead = results->mappedBlocks;
    hal_block_t *prevBlock = NULL;
    while (blockHead != NULL) {
        hal_block_t *next = blockHead->next;
        if (strcmp(blockHead->qChrom, qChrom) != 0) {
            free(blockHead->qSequence);
            free(blockHead->tSequence);
            free(blockHead->qChrom);
            free(blockHead);
        } else {
            if (prevBlock != NULL) {
                prevBlock->next = blockHead;
            } else {
                results->mappedBlocks = blockHead;
            }
            prevBlock = blockHead;
        };
        blockHead = next;
    }
    if (prevBlock != NULL) {
        prevBlock->next = NULL;
    } else {
        results->mappedBlocks = NULL;
    }

    // Filter out any target dupes not corresponding to qChrom.
    hal_target_dupe_list_t *targetDupeHead = results->targetDupeBlocks;
    hal_target_dupe_list_t *prevTargetDupe = NULL;
    while (targetDupeHead != NULL) {
        hal_target_dupe_list_t *next = targetDupeHead->next;
        if (strcmp(targetDupeHead->qChrom, qChrom) != 0) {
            while (targetDupeHead->tRange != NULL)
            {
                hal_target_range_t* nextRange = targetDupeHead->tRange->next;
                free(targetDupeHead->tRange);
                targetDupeHead->tRange = nextRange;
            }
            free(targetDupeHead->qChrom);
            free(targetDupeHead);
        } else {
            if (prevTargetDupe != NULL) {
                prevTargetDupe->next = targetDupeHead;
            } else {
                results->targetDupeBlocks = targetDupeHead;
            }
            prevTargetDupe = targetDupeHead;
        }
        targetDupeHead = next;
    }
    if (prevTargetDupe != NULL) {
        prevTargetDupe->next = NULL;
    } else {
        results->targetDupeBlocks = NULL;
    }
    return results;
}

extern "C" hal_int_t halGetMAF(FILE* outFile,
                               int halHandle, 
                               hal_species_t* qSpeciesNames,
                               char* tSpecies,
                               char* tChrom,
                               hal_int_t tStart, 
                               hal_int_t tEnd,
                               int maxRefGap,
                               int maxBlockLength,
                               int doDupes,
                               char **errStr)
{
  halLock();
  hal_int_t numBytes = 0;
  try
  {
    hal_int_t rangeLength = tEnd - tStart;
    if (rangeLength < 0)
    {
        halUnlock();
        handleError("halGetMAF invalid query range [" + std::to_string(tStart) + "," + std::to_string(tEnd) + ")", errStr);
        return -1;
    }
    const Alignment* alignment = getExistingAlignment(halHandle, hal_size_t(0), true);

    set<const Genome*> qGenomeSet;
    for (hal_species_t* qSpecies = qSpeciesNames; qSpecies != NULL;
         qSpecies = qSpecies->next)
    {
        checkGenomes(halHandle, alignment, qSpecies->name, tSpecies, tChrom);
      const Genome* qGenome = alignment->openGenome(qSpecies->name);
      qGenomeSet.insert(qGenome);
    }
    
    const Genome* tGenome = alignment->openGenome(tSpecies);
    const Sequence* tSequence = tGenome->getSequence(tChrom);

    hal_index_t myEnd = tEnd > 0 ? tEnd : tSequence->getSequenceLength();
    hal_index_t absStart = tSequence->getStartPosition() + tStart;
    hal_index_t absEnd = tSequence->getStartPosition() + myEnd - 1;
    if (absStart > absEnd)
    {
        halUnlock();
        handleError("halGetMAF invalid range", errStr);
        return -1;
    }
    if (absEnd > tSequence->getEndPosition())
    {
        halUnlock();
        handleError("halGetMAF target end position outside of target sequence", errStr);
        return -1;
    }

    stringstream mafBuffer;
    MafExport mafExport;
    mafExport.setNoDupes(doDupes == 0);
    mafExport.setUcscNames(true);
    mafExport.setMaxRefGap(hal_size_t(maxRefGap));
    mafExport.setMaxBlockLength(hal_index_t(maxBlockLength));

    mafExport.convertSegmentedSequence(mafBuffer, alignment, tGenome, 
                                       absStart, 1 + absEnd - absStart,
                                       qGenomeSet);
    // if these buffers are very big, it is my intuition that 
    // chopping them up and, saw, only converting and writing a few kb
    // at a time would be more efficient... but i think that's the least
    // of our problems right now.  
    string mafStringBuffer = mafBuffer.str();
    if (mafStringBuffer.empty() == false)
    {
      numBytes = (hal_int_t)fwrite(mafStringBuffer.c_str(), 
                                   mafStringBuffer.length(), 
                                   sizeof(char), outFile);
    }                 
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetMAF error writing MAF blocks: " + string(e.what()), errStr);
    return -1;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetMAF error writing MAF blocks: unknown exception", errStr);
    return -1;
  }
  halUnlock();
  return numBytes;
}

extern "C" struct hal_species_t *halGetSpecies(int halHandle, char **errStr)
{
  halLock();
  hal_species_t* head = NULL;
  try
  {
    // read the lowest level of detail because it's fastest
      const Alignment* alignment = getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
                                                       false);
    hal_species_t* prev = NULL;
    if (alignment->getNumGenomes() > 0)
    {
      string rootName = alignment->getRootName();
      deque<string> bfQueue(1, rootName);
      while (!bfQueue.empty())
      {
        string name = bfQueue.back();
        bfQueue.pop_back();
        const Genome* genome = alignment->openGenome(name);
        
        hal_species_t* cur = (hal_species_t*)calloc(1, sizeof(hal_species_t));
        cur->next = NULL;
        cur->name = copyCString(name);
        cur->length = genome->getSequenceLength();
        cur->numChroms = genome->getNumSequences();
        if (name == rootName)
        {
          cur->parentName = NULL;
          cur->parentBranchLength = 0;
        }
        else
        {
          cur->parentName = copyCString(alignment->getParentName(name));
          cur->parentBranchLength = 
             alignment->getBranchLength(cur->parentName, name);
        }
        if (head == NULL)
        {
          head = cur;
        }
        else
        {
          prev->next = cur;
        }
        prev = cur;
        
        vector<string> childNames = alignment->getChildNames(name);
        for (size_t i = 0; i < childNames.size(); ++i)
        {
          bfQueue.push_front(childNames[i]);
        }
      }      
    }
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetSpecies: " + string(e.what()), errStr);
    return NULL;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetSpecies: unknown exception", errStr);
    return NULL;
  }
  halUnlock();
  return head;
}

extern "C" struct hal_species_t *halGetPossibleCoalescenceLimits(int halHandle,
                                                                 const char *qSpecies,
                                                                 const char *tSpecies,
                                                                 char **errStr) {
  halLock();
  hal_species_t* head = NULL;
  try
  {
    // read the lowest level of detail because it's fastest
      const Alignment* alignment = getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
                                                        false);
    hal_species_t* prev = NULL;
    const Genome *qGenome = alignment->openGenome(qSpecies);
    const Genome *tGenome = alignment->openGenome(tSpecies);
    set<const Genome*> inputSet;
    inputSet.insert(qGenome);
    inputSet.insert(tGenome);
    const Genome *mrca = getLowestCommonAncestor(inputSet);

    // Get a list of all ancestors of the MRCA.
    const Genome *curGenome = mrca;
    do {
      hal_species_t* cur = (hal_species_t*)calloc(1, sizeof(hal_species_t));
      cur->next = NULL;
      cur->name = copyCString(curGenome->getName());
      cur->length = curGenome->getSequenceLength();
      cur->numChroms = curGenome->getNumSequences();
      if (curGenome->getName() == alignment->getRootName())
      {
        cur->parentName = NULL;
        cur->parentBranchLength = 0;
      }
      else
      {
        cur->parentName = copyCString(curGenome->getParent()->getName());
        cur->parentBranchLength = 
          alignment->getBranchLength(cur->parentName, cur->name);
      }
      if (head == NULL)
      {
        head = cur;
      }
      else
      {
        prev->next = cur;
      }
      prev = cur;
    } while ((curGenome = curGenome->getParent()) != NULL);
  }
  catch(exception& e)
  {
      halUnlock();
      handleError("halGetPossibleCoalescenceLimits: " + string(e.what()), errStr);
      return NULL;
  }
  catch(...)
  {
      halUnlock();
      handleError("halGetPossibleCoalescenceLimits: unknown exception", errStr);
      return NULL;
  }
  halUnlock();
  return head;
}

extern "C" void halFreeSpeciesList(struct hal_species_t *species)
{
  while (species != NULL) {
    free(species->name);
    free(species->parentName);
    struct hal_species_t *nextSpecies = species->next;
    free(species);
    species = nextSpecies;
  }
}

extern "C" struct hal_chromosome_t *halGetChroms(int halHandle,
                                                 char* speciesName,
                                                 char **errStr)
{
  halLock();
  hal_chromosome_t* head = NULL;
  try
  {
    // read the lowest level of detail because it's fastest
      const Alignment* alignment = getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
                                                        false);

    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
    halUnlock();
    handleError("halGetChroms: species with name " + string(speciesName)
                + " not found in alignment with handle "
                + std::to_string(halHandle), errStr);
    return NULL;
    }

    hal_chromosome_t* prev = NULL;
    if (genome->getNumSequences() > 0)
    {
      SequenceIteratorPtr seqIt = genome->getSequenceIterator();
      for (; not seqIt->atEnd(); seqIt->toNext())
      {
        const Sequence* sequence = seqIt->getSequence();
        
        hal_chromosome_t* cur = 
           (hal_chromosome_t*)calloc(1, sizeof(hal_chromosome_t));
        cur->next = NULL;
        cur->name = copyCString(sequence->getName());
        cur->length = sequence->getSequenceLength();

        if (head == NULL)
        {
          head = cur;
        }
        else
        {
          prev->next = cur;
        }
        prev = cur;
      }
    }
  }
  catch(exception& e)
  {
      halUnlock();
      handleError("halGetChroms: " + string(e.what()) , errStr);
      return NULL;
  }
  catch(...)
  {
      halUnlock();
      handleError("halGetChroms: unknown exception", errStr);
      return NULL;
  }
  halUnlock();
  return head;
}

extern "C" char *halGetDna(int halHandle,
                           char* speciesName, char* chromName, 
                           hal_int_t start, hal_int_t end,
                           char **errStr)
{
  halLock();
  char* dna = NULL;
  try
  {
      const Alignment* alignment = getExistingAlignment(halHandle, 0, true);
    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
        throw hal_exception("halGetDna: species with name " + string(speciesName) + " not found in alignment with handle "
                            + std::to_string(halHandle));
    }
    const Sequence* sequence = genome->getSequence(chromName);
    if (sequence == NULL)
        {
            halUnlock();
            handleError("halGetDna: chromosome with name " + string(chromName) + " not found in species "
                        + speciesName, errStr);
            return NULL;
    }    
    if (start > end || end > (hal_index_t)sequence->getSequenceLength())
    {
    halUnlock();
    handleError("halGetDna: specified range [" + std::to_string(start) + "," + std::to_string(end) + ") is invalid "
                + "for chromsome " + chromName + " in species " + speciesName
                + " which is of length " + std::to_string(sequence->getSequenceLength()), errStr);
    return NULL;
    }
    
    string buffer;
    sequence->getSubString(buffer, start, end - start);
    dna = copyCString(buffer);
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetDna: " + string(e.what()), errStr);
    return NULL;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetDna: unknown exception", errStr);
    return NULL;
  }
  halUnlock();
  return dna;
}

extern "C" hal_int_t halGetMaxLODQueryLength(int halHandle, char **errStr)
{
  halLock();
  hal_int_t ret = 0;
  try
  {
    HandleMap::iterator mapIt = handleMap.find(halHandle);
    if (mapIt == handleMap.end())
    {
    halUnlock();
    handleError("halGetMaxLODQueryLength error getting Max LOD Query Length.  handle " 
                + std::to_string(halHandle) + ": not found", errStr);
    return -1;
    }
    ret = (hal_int_t)mapIt->second.second->getMaxQueryLength();
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetMaxLODQueryLength: " + string(e.what()), errStr);
    return -1;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetMaxLODQueryLength: unknown exception", errStr);
    return -1;
  }
  halUnlock();
  return ret;
}

static void checkHandle(int handle)
{
  HandleMap::iterator mapIt = handleMap.find(handle);
  if (mapIt == handleMap.end())
  {
      throw hal_exception("Handle " + std::to_string(handle) + "not found in alignment map");
  }
  if (mapIt->second.second.get() == NULL)
  {
      throw hal_exception( "Handle " + std::to_string(handle) + "points to NULL alignment");
  }
}

static void checkGenomes(int halHandle, 
                         const Alignment* alignment, const string& qSpecies,
                         const string& tSpecies, const string& tChrom)
{
  const Genome* qGenome = alignment->openGenome(qSpecies);
  if (qGenome == NULL)
  {
      throw hal_exception("Query species " + qSpecies + " not found in alignment with "
                          + "handle " + std::to_string(halHandle));
  }
  const Genome* tGenome = alignment->openGenome(tSpecies);
  if (tGenome == NULL)
  {
    throw hal_exception("Reference species " + tSpecies + " not found in alignment with "
                        + "handle " + std::to_string(halHandle));
  }

  const Sequence* tSequence = tGenome->getSequence(tChrom);
  if (tSequence == NULL)
  {
      throw hal_exception("Unable to locate sequence " + tChrom + " in genome " + tSpecies);
  }
}


static const Alignment* getExistingAlignment(int handle, hal_size_t queryLength,
                                              bool needDNASequence)
{
  checkHandle(handle);
  HandleMap::iterator mapIt = handleMap.find(handle);
  return mapIt->second.second->getAlignment(queryLength, needDNASequence);
}

static bool isAlignmentLod0(int handle, hal_size_t queryLength)
{
  checkHandle(handle);
  HandleMap::iterator mapIt = handleMap.find(handle);
  return mapIt->second.second->isLod0(queryLength);
}

static char* copyCString(const string& inString)
{
  char* outString = (char*)malloc(inString.length() + 1);
  strcpy(outString, inString.c_str());
  return outString;
}

static hal_block_results_t* readBlocks(const Alignment* seqAlignment,
                                       const Sequence* tSequence,
                                       hal_index_t absStart, hal_index_t absEnd,
                                       bool tReversed,
                                       const Genome* qGenome, bool getSequenceString,
                                       bool doDupes, bool doTargetDupes, 
                                       bool doAdjes, const char *coalescenceLimitName)
{
  const Genome* tGenome = tSequence->getGenome();
  string qGenomeName = qGenome->getName();
  hal_block_t* prev = NULL;
  BlockMapper blockMapper;
  if (qGenome == tGenome && coalescenceLimitName == NULL)
  {
    // By default, for self-alignment tracks, walk all the way back to
    // the root finding paralogies.
    const Alignment *alignment = tGenome->getAlignment();
    const Genome *root = alignment->openGenome(alignment->getRootName());
    blockMapper.init(tGenome, qGenome, absStart, absEnd,
                     tReversed, doDupes, 0, doAdjes, root);
  }
  else
  {
    // Just use the standard parameters, with coalescence limit set to
    // the genome represented by coalescenceLimitName, or the MRCA if
    // NULL.
    const Alignment *alignment = tGenome->getAlignment();
    const Genome *coalescenceLimit = NULL;

    if (coalescenceLimitName != NULL) {
      coalescenceLimit = alignment->openGenome(coalescenceLimitName);
      if (coalescenceLimit == NULL) {
        throw hal_exception("Could not find coalescence limit "
                            + string(coalescenceLimitName) + " in alignment");
      }
    }

    blockMapper.init(tGenome, qGenome, absStart, absEnd,
                     tReversed, doDupes, 0, doAdjes, coalescenceLimit);
  }
  blockMapper.map();
  MappedSegmentSet paraSet;
  hal_size_t totalLength = 0;
  hal_size_t reversedLength = 0;
  if (doDupes == true && qGenome != tGenome)
  {
    blockMapper.extractReferenceParalogies(paraSet);
  }
  MappedSegmentSet& segMap = blockMapper.getMap();
  vector<MappedSegmentPtr> fragments;
  set<hal_index_t> queryCutSet;
  set<hal_index_t> targetCutSet;
  targetCutSet.insert(blockMapper.getAbsRefFirst());
  targetCutSet.insert(blockMapper.getAbsRefLast());

  hal_block_results_t* results = 
     (hal_block_results_t*)calloc(1, sizeof(hal_block_results_t));

  for (MappedSegmentSet::iterator segMapIt = segMap.begin();
       segMapIt != segMap.end(); ++segMapIt)
  {
    assert((*segMapIt)->getSource()->getReversed() == false);
    hal_block_t* cur = (hal_block_t*)calloc(1, sizeof(hal_block_t));
    if (results->mappedBlocks == NULL)
    {
      results->mappedBlocks = cur;
    }
    else
    {
      prev->next = cur;
    }
    BlockMapper::extractSegment(segMapIt, paraSet, fragments, &segMap,
                                targetCutSet, queryCutSet);
    readBlock(seqAlignment, cur, fragments, getSequenceString, qGenomeName);
    totalLength += cur->size;
    reversedLength += cur->strand == '-' ? cur->size : 0;
    prev = cur;
  }
  if (!paraSet.empty())
  {
    results->targetDupeBlocks = processTargetDupes(blockMapper, paraSet);
  }
  if (doTargetDupes == false && results->targetDupeBlocks != NULL)
  {
    halFreeTargetDupeLists(results->targetDupeBlocks);
    results->targetDupeBlocks = NULL;
  }
  return results;
}

static void readBlock(const Alignment* seqAlignment,
                      hal_block_t* cur,  
                      vector<MappedSegmentPtr>& fragments, 
               bool getSequenceString, const string& genomeName)
{
  MappedSegmentPtr firstQuerySeg = fragments.front();
  MappedSegmentPtr lastQuerySeg = fragments.back();
  const SlicedSegment* firstRefSeg = firstQuerySeg->getSource();
  const SlicedSegment* lastRefSeg = lastQuerySeg->getSource();
  const Sequence* qSequence = firstQuerySeg->getSequence();
  const Sequence* tSequence = firstRefSeg->getSequence(); 
  assert(qSequence == lastQuerySeg->getSequence());
  assert(tSequence == lastRefSeg->getSequence());
  assert(firstRefSeg->getReversed() == false);
  assert(lastRefSeg->getReversed() == false);
  
  cur->next = NULL;

  string seqBuffer = qSequence->getName();
  string qDnaBuffer;
  string tDnaBuffer;
  size_t prefix = 
     seqBuffer.find(genomeName + '.') != 0 ? 0 : genomeName.length() + 1;
  cur->qChrom = (char*)malloc(seqBuffer.length() + 1 - prefix);
  strcpy(cur->qChrom, seqBuffer.c_str() + prefix);

  cur->tStart = std::min(std::min(firstRefSeg->getStartPosition(), 
                                  firstRefSeg->getEndPosition()),
                         std::min(lastRefSeg->getStartPosition(),
                                  lastRefSeg->getEndPosition()));
  cur->tStart -= tSequence->getStartPosition();

  cur->qStart = std::min(std::min(firstQuerySeg->getStartPosition(), 
                                  firstQuerySeg->getEndPosition()),
                         std::min(lastQuerySeg->getStartPosition(),
                                  lastQuerySeg->getEndPosition()));
  cur->qStart -= qSequence->getStartPosition();

  hal_index_t tEnd = std::max(std::max(firstRefSeg->getStartPosition(), 
                                       firstRefSeg->getEndPosition()),
                              std::max(lastRefSeg->getStartPosition(),
                                       lastRefSeg->getEndPosition()));
  tEnd -= tSequence->getStartPosition();

  assert(cur->tStart >= 0);
  assert(cur->qStart >= 0);

  assert(firstRefSeg->getLength() == firstQuerySeg->getLength());
  cur->size = 1 + tEnd - cur->tStart;
  cur->strand = firstQuerySeg->getReversed() ? '-' : '+';
  cur->tSequence = NULL;
  cur->qSequence = NULL;
  if (getSequenceString != 0)
  {
    const Genome* qSeqGenome = 
       seqAlignment->openGenome(qSequence->getGenome()->getName());
    if (qSeqGenome == NULL)
    {
      throw hal_exception("Unable to open genome " + qSequence->getGenome()->getName() 
                          + " for DNA sequence extraction");
    }
    const Sequence* qSeqSequence = qSeqGenome->getSequence(qSequence->getName());
    if (qSeqSequence == NULL)
    {
      throw hal_exception("Unable to open sequence " + qSequence->getName() 
                          + " for DNA sequence extraction");
    }

    const Genome* tSeqGenome = 
       seqAlignment->openGenome(tSequence->getGenome()->getName());
    if (tSeqGenome == NULL)
    {
      throw hal_exception("Unable to open genome " + tSequence->getGenome()->getName() 
                          + " for DNA sequence extraction");
    }

    const Sequence* tSeqSequence = tSeqGenome->getSequence(tSequence->getName());
    if (tSeqSequence == NULL)
    {
      throw hal_exception("Unable to open sequence " + tSequence->getName() 
                          + " for DNA sequence extraction");
    }
    
    qSeqSequence->getSubString(qDnaBuffer, cur->qStart, cur->size);
    tSeqSequence->getSubString(tDnaBuffer, cur->tStart, cur->size);
    if (cur->strand == '-')
    {
      reverseComplement(qDnaBuffer);
    }
    cur->qSequence = (char*)malloc(qDnaBuffer.length() * sizeof(char) + 1);
    cur->tSequence = (char*)malloc(tDnaBuffer.length() * sizeof(char) + 1);
    strcpy(cur->qSequence, qDnaBuffer.c_str());
    strcpy(cur->tSequence, tDnaBuffer.c_str());
  }
}

struct CStringLess { bool operator()(const char* s1, const char* s2) const {
  return strcmp(s1, s2) < 0;}
};

struct DupeIdLess { bool operator()(const hal_target_dupe_list_t* d1,
                                    const hal_target_dupe_list_t* d2) const {
  if (d1 != NULL && d2 != NULL)
  {
    if (d1->id == d2->id)
    {
      // hack to keep ties in reverse sorted order based on tStart
      // this is so that the prepend in processTargetDupes when the 
      // lists gets merged keeps the tStarts in ascending order
      return d1->tRange->tStart > d2->tRange->tStart;
    }
    return d1->id < d2->id;
  }
  return d1 != NULL && d2 == NULL;
}};

struct DupeStartLess { bool operator()(const hal_target_dupe_list_t* d1,
                                       const hal_target_dupe_list_t* d2) const {
  if (d1 != NULL && d2 != NULL)
  {
    return d1->tRange->tStart < d2->tRange->tStart;
  }
  return d1 != NULL && d2 == NULL;
}};

static hal_target_dupe_list_t* processTargetDupes(BlockMapper& blockMapper,
                                                  MappedSegmentSet& paraSet)
{
  vector<hal_target_dupe_list_t*> tempList;
  vector<MappedSegmentPtr> fragments;
  MappedSegmentSet emptySet;
  set<hal_index_t> queryCutSet;
  set<hal_index_t> targetCutSet;
  targetCutSet.insert(blockMapper.getAbsRefFirst());
  targetCutSet.insert(blockMapper.getAbsRefLast());

  // make a dupe list for each merged segment
  for (MappedSegmentSet::iterator segMapIt = paraSet.begin();
       segMapIt != paraSet.end(); ++segMapIt)
  {
    assert((*segMapIt)->getSource()->getReversed() == false);
    hal_target_dupe_list_t* cur = (hal_target_dupe_list_t*)calloc(
      1, sizeof(hal_target_dupe_list_t));
    BlockMapper::extractSegment(segMapIt, emptySet, fragments, &paraSet,
                                targetCutSet, queryCutSet);
    //fragments.clear(); fragments.push_back(*segMapIt);
    readTargetRange(cur, fragments);
    tempList.push_back(cur);
  }
  
  cleanTargetDupesList(tempList);

  map<const char*, hal_int_t, CStringLess> idMap;
  pair<map<const char*, hal_int_t, CStringLess>::iterator, bool> idRes;
  vector<hal_target_dupe_list_t*>::iterator i;
  vector<hal_target_dupe_list_t*>::iterator j;
  vector<hal_target_dupe_list_t*>::iterator k;
  for (i = tempList.begin(); i != tempList.end(); i = j)
  {
    j = i;
    ++j;
    // merge dupe lists with overlapping ids by prepending j's
    // target range to i
    // (recall we've stuck query start coordinates into the id)
    while (j != tempList.end() && (*i)->id == (*j)->id)
    {
      assert((*i)->next == NULL);
      assert((*j)->next == NULL);
      assert(j != i);
      k = j;
      ++k;
      assert((*j)->tRange->size == (*i)->tRange->size);
      (*j)->tRange->next = (*i)->tRange;
      (*i)->tRange = (*j)->tRange;
      (*j)->tRange = NULL;
      // DupeIdLess sorting comparator should ensure that these 
      // elemens are added in the proper order so that trange list
      // stays sorted here
      assert((*i)->tRange->tStart <= (*i)->tRange->next->tStart);
      halFreeTargetDupeLists(*j);
      *j = NULL;
      j = k;
    }

    idRes = idMap.insert(pair<const char*, hal_int_t>((*i)->qChrom, 0));
    (*i)->id = idRes.first->second++;
  }
 
  std::sort(tempList.begin(), tempList.end(), DupeStartLess());
  mergeCompatibleDupes(tempList);

  hal_target_dupe_list_t* head = NULL;
  hal_target_dupe_list_t* prev = NULL;
  for (i = tempList.begin(); i != tempList.end() && *i != NULL; ++i)
  {
    if (head == NULL)
    {
      head = *i;
    }
    else 
    {
      prev->next = *i;
    }
    prev = *i;
  }
  
  return head;
}

static void readTargetRange(hal_target_dupe_list_t* cur,
                            vector<MappedSegmentPtr>& fragments)
{
  MappedSegmentPtr firstQuerySeg = fragments.front();
  MappedSegmentPtr lastQuerySeg = fragments.back();
  const SlicedSegment* firstRefSeg = firstQuerySeg->getSource();
  const SlicedSegment* lastRefSeg = lastQuerySeg->getSource();
  const Sequence* qSequence = firstQuerySeg->getSequence();
  const Sequence* tSequence = firstRefSeg->getSequence(); 
  assert(qSequence == lastQuerySeg->getSequence());
  assert(tSequence == lastRefSeg->getSequence());
  assert(firstRefSeg->getReversed() == false);
  assert(lastRefSeg->getReversed() == false);

  cur->tRange = (hal_target_range_t*)calloc(1, sizeof(hal_target_range_t));
 
  cur->tRange->tStart = std::min(std::min(firstRefSeg->getStartPosition(), 
                                          firstRefSeg->getEndPosition()),
                                 std::min(lastRefSeg->getStartPosition(),
                                          lastRefSeg->getEndPosition()));
  cur->tRange->tStart -= tSequence->getStartPosition();

  hal_index_t qStart = std::min(std::min(firstQuerySeg->getStartPosition(), 
                                         firstQuerySeg->getEndPosition()),
                                std::min(lastQuerySeg->getStartPosition(),
                                         lastQuerySeg->getEndPosition()));
  qStart -= qSequence->getStartPosition();

  hal_index_t tEnd = std::max(std::max(firstRefSeg->getStartPosition(), 
                                       firstRefSeg->getEndPosition()),
                              std::max(lastRefSeg->getStartPosition(),
                                       lastRefSeg->getEndPosition()));
  tEnd -= tSequence->getStartPosition();
  
  cur->tRange->size = 1 + tEnd - cur->tRange->tStart;

  // use query start as proxy for unique id right now
  cur->id = qStart;

  string qSeqName = qSequence->getName();
  cur->qChrom = (char*)malloc(qSeqName.length() * sizeof(char) + 1);
  strcpy(cur->qChrom, qSeqName.c_str());
}

static void cleanTargetDupesList(vector<hal_target_dupe_list_t*>& dupeList)
{
  // sort based on target start coordinate
  std::sort(dupeList.begin(), dupeList.end(), DupeStartLess());
  
  map<hal_int_t, hal_int_t> idMap;
  map<hal_int_t, hal_int_t>::iterator idMapIt;

  size_t deadCount = 0;

  vector<hal_target_dupe_list_t*>::iterator i;
  vector<hal_target_dupe_list_t*>::iterator j;
  vector<hal_target_dupe_list_t*>::iterator k;
  for (j = dupeList.begin(); j != dupeList.end(); j = k)
  {
    k = j;
    ++k;
    if (j == dupeList.begin() || 
        // note that we should eventually clean up overlapping
        // segments but don't have the energy to do right now, and
        // have yet to see it actually happen in an example
        strcmp((*j)->qChrom, (*i)->qChrom) != 0 ||
        (*j)->tRange->tStart != (*i)->tRange->tStart ||
        (*j)->tRange->size != (*i)->tRange->size)
    {
      i = j;
    }
    else
    {
      idMap.insert(pair<hal_int_t, hal_int_t>((*j)->id, (*i)->id));
      // insert this one to make sure we dont remap it down the road
      idMap.insert(pair<hal_int_t, hal_int_t>((*i)->id, (*i)->id));
      halFreeTargetDupeLists(*j);
      *j = NULL;
      ++deadCount;
    }
  }

  for (i = dupeList.begin(); i != dupeList.end(); ++i)
  {
    if (*i != NULL)
    {
      idMapIt = idMap.find((*i)->id);
      if (idMapIt != idMap.end())
      {
        (*i)->id = idMapIt->second;
      }
    }
  }

  // sort based on query coordinate
  std::sort(dupeList.begin(), dupeList.end(), DupeIdLess());

  // sort puts NULL items at end. we resize them out.
  dupeList.resize(dupeList.size() - deadCount);
}

static bool areRangesCompatible(hal_target_range_t* range1,
                                hal_target_range_t* range2)
{
  for (; range1 != NULL; range1 = range1->next, range2 = range2->next)
  {
    if ((range1->next == NULL) != (range2->next == NULL) ||
        range2->tStart != range1->tStart + range1->size ||
        (range1->next != NULL && 
         range1->next->tStart - (range1->tStart + range1->size) < 
         range2->size))
    {
      return false;
    }        
  }
  return true;
}

static void mergeCompatibleDupes(vector<hal_target_dupe_list_t*>& dupeList)
{
  vector<hal_target_dupe_list_t*>::iterator i;
  vector<hal_target_dupe_list_t*>::iterator j;
  vector<hal_target_dupe_list_t*>::iterator k;
  
  for (i = dupeList.begin(); i != dupeList.end() && *i != NULL; i = j)
  {
    j = i;
    bool compatible = true;
    for (++j; j != dupeList.end() && *j != NULL && compatible; ++j)
    {
      hal_target_range_t* range1 = (*i)->tRange;
      hal_target_range_t* range2 = (*j)->tRange;
      compatible = areRangesCompatible(range1, range2);
      if (compatible == true)
      {
        range2 = (*j)->tRange;
        for (range1 = (*i)->tRange; range1 != NULL; range1 = range1->next) 
        {
          assert(range2 != NULL);
          range1->size += range2->size;
          range2 = range2->next;
        } 
        halFreeTargetDupeLists(*j);
        *j = NULL;
      }
    }
  }
}

extern "C" struct hal_metadata_t *halGetGenomeMetadata(int halHandle,
                                                       const char *genomeName,
                                                       char **errStr)
{
  halLock();
  struct hal_metadata_t *ret = NULL;
  try {
      const Alignment* alignment = getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
                                                        false);

    const Genome *genome = alignment->openGenome(genomeName);
    if (genome == NULL) {
        throw hal_exception("Genome " + string(genomeName) + " not found in alignment");
    }
    const MetaData *metaData = genome->getMetaData();
    const map<string, string> metaDataMap = metaData->getMap();

    hal_metadata_t *prevMetadata = NULL;
    for (map<string, string>::const_iterator i = metaDataMap.begin();
         i != metaDataMap.end(); i++) {
      hal_metadata_t *curMetadata = (hal_metadata_t *) calloc(1, sizeof(hal_metadata_t));
      curMetadata->key = copyCString(i->first.c_str());
      curMetadata->value = copyCString(i->second.c_str());
      if (prevMetadata != NULL) {
        prevMetadata->next = curMetadata;
      } else {
        ret = curMetadata;
      }
      prevMetadata = curMetadata;
    }
  }
  catch(exception& e)
  {
    halUnlock();
    handleError("halGetGenomeMetadata: " + string(e.what()), errStr);
    return NULL;
  }
  catch(...)
  {
    halUnlock();
    handleError("halGetGenomeMetadata: unknown exception", errStr);
    return NULL;
  }
  halUnlock();
  return ret;
}

extern "C" void halFreeMetadataList(struct hal_metadata_t *metadata) {
  while (metadata != NULL) {
    struct hal_metadata_t *next = metadata->next;
    free(metadata->key);
    free(metadata->value);
    free(metadata);
    metadata = next;
  }
}

extern "C" void halFreeChromList(struct hal_chromosome_t *chroms) {
  while (chroms != NULL) {
    free(chroms->name);
    struct hal_chromosome_t *tmp = chroms->next;
    free(chroms);
    chroms = tmp;
  }
}
