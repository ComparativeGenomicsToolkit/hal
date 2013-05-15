/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <map>
#include <sstream>
#include <limits>
#include <stdlib.h>
#include <string.h>
#include "hal.h"
#include "halChain.h"
#include "halBlockViz.h"
#include "halBlockMapper.h"
#include "halLodManager.h"

#ifdef ENABLE_UDC
#include <pthread.h>
static pthread_mutex_t HAL_MUTEX;
#define HAL_LOCK pthread_mutex_lock(&HAL_MUTEX);
#define HAL_UNLOCK pthread_mutex_unlock(&HAL_MUTEX);
#else
#define HAL_LOCK
#define HAL_UNLOCK
#endif

using namespace std;
using namespace hal;

typedef map<int, pair<string, LodManagerPtr> > HandleMap;
static HandleMap handleMap;

static int halOpenLodOrHal(char* inputPath, bool isLod);
static void checkHandle(int handle);
static void checkGenomes(int halHandle, 
                         AlignmentConstPtr alignment, const string& qSpecies,
                         const string& tSpecies, const string& tChrom);

static AlignmentConstPtr getExistingAlignment(int handle,
                                              hal_size_t queryLength,
                                              bool needSequence);
static char* copyCString(const string& inString);

static hal_block_results_t* readBlocks(const Sequence* tSequence,
                                       hal_index_t absStart, hal_index_t absEnd,
                                       const Genome* qGenome, 
                                       bool getSequenceString,
                                       bool doDupes);

static void readBlock(hal_block_t* cur, 
                      vector<MappedSegmentConstPtr>& fragments,                                        bool getSequenceString, const string& genomeName);

static hal_target_dupe_list_t* processTargetDupes(BlockMapper& blockMapper,
                                                  BlockMapper::MSSet& paraSet);

static void readTargetRange(hal_target_dupe_list_t* cur,
                            vector<MappedSegmentConstPtr>& fragments);


extern "C" int halOpenLOD(char *lodFilePath)
{
  bool isHal = lodFilePath && strlen(lodFilePath) > 4 &&
     strcmp(lodFilePath + strlen(lodFilePath) - 4, ".hal") == 0;

  return halOpenLodOrHal(lodFilePath, !isHal);
}

extern "C" int halOpen(char* halFilePath)
{
  return halOpenLodOrHal(halFilePath, false);
}

int halOpenLodOrHal(char* inputPath, bool isLod)
{
  HAL_LOCK
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
    cerr << "Exception caught: " << e.what() << endl;
    handle = -1;
  }
  catch(...)
  {
    cerr << "Error opening " << inputPath << endl;
    handle = -1;
  }
  HAL_UNLOCK
  return handle;
}

extern "C" int halClose(int handle)
{
  HAL_LOCK
  int ret = 0;
  try
  {
    HandleMap::iterator mapIt = handleMap.find(handle);
    if (mapIt == handleMap.end())
    {
      stringstream ss;
      ss << "error closing handle " << handle << ": not found";
      throw hal_exception(ss.str());
    }
    handleMap.erase(mapIt);
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    ret = -1;
  }
  catch(...)
  {
    cerr << "Error closing " << handle << endl;
    ret = -1;
  }
  HAL_UNLOCK
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
    free(head->sequence);
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
                                                      int getSequenceString,
                                                      int doDupes)
{
  HAL_LOCK
  hal_block_results_t* results = NULL;
  try
  {
    hal_int_t rangeLength = tEnd - tStart;
    if (rangeLength < 0)
    {
      stringstream ss;
      ss << "Invalid query range [" << tStart << "," << tEnd << ").";
      throw hal_exception(ss.str());
    }
    AlignmentConstPtr alignment = 
       getExistingAlignment(halHandle, hal_size_t(rangeLength), 
                            getSequenceString != 0);
    checkGenomes(halHandle, alignment, qSpecies, tSpecies, tChrom);

    const Genome* qGenome = alignment->openGenome(qSpecies);
    const Genome* tGenome = alignment->openGenome(tSpecies);
    const Sequence* tSequence = tGenome->getSequence(tChrom);
    // cactus pipeline presently adds species name as prefix of 
    // sequence name.  check if this caused confusion
    string sequenceName;
    if (tSequence == NULL)
    {
      sequenceName = tGenome->getName();
      sequenceName += '.';
      sequenceName += tChrom;
      tSequence = tGenome->getSequence(sequenceName);
    }

    hal_index_t myEnd = tEnd > 0 ? tEnd : tSequence->getSequenceLength();
    hal_index_t absStart = tSequence->getStartPosition() + tStart;
    hal_index_t absEnd = tSequence->getStartPosition() + myEnd - 1;
    if (absStart > absEnd)
    {
      throw hal_exception("Invalid range");
    }
    // We now know the query length so we can do a proper lod query
    if (tEnd == 0)
    {
      alignment = getExistingAlignment(halHandle, absEnd - absStart, 
                                       getSequenceString != 0);
      checkGenomes(halHandle, alignment, qSpecies, tSpecies, tChrom);
      qGenome = alignment->openGenome(qSpecies);
      tGenome = alignment->openGenome(tSpecies);
      tSequence = tGenome->getSequence(tSequence->getName());
    }

    results = readBlocks(tSequence, absStart, absEnd, qGenome,
                         getSequenceString != 0, doDupes != 0);
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    results = NULL;
  }
  catch(...)
  {
    cerr << "Error in hal block query";
    results = NULL;
  }
  HAL_UNLOCK
  return results;
}

extern "C" struct hal_species_t *halGetSpecies(int halHandle)
{
  HAL_LOCK
  hal_species_t* head = NULL;
  try
  {
    // read the lowest level of detail because it's fastest
    AlignmentConstPtr alignment = 
       getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
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
        
        alignment->closeGenome(genome);
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
    cerr << "Exception caught: " << e.what() << endl;
    head = NULL;
  }
  catch(...)
  {
    cerr << "Error in hal get species";
    head = NULL;
  }
  HAL_UNLOCK
  return head;
}

extern "C" struct hal_chromosome_t *halGetChroms(int halHandle,
                                                 char* speciesName)
{
  HAL_LOCK
  hal_chromosome_t* head = NULL;
  try
  {
    // read the lowest level of detail because it's fastest
    AlignmentConstPtr alignment = 
       getExistingAlignment(halHandle, numeric_limits<hal_size_t>::max(), 
                            false);

    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
      stringstream ss;
      ss << "Species with name " << speciesName << " not found in alignment "
         << "with handle " << halHandle;
      throw hal_exception(ss.str());
    }

    hal_chromosome_t* prev = NULL;
    if (genome->getNumSequences() > 0)
    {
      SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
      SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
      for (; seqIt != seqEnd; seqIt->toNext())
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
    cerr << "Exception caught: " << e.what() << endl;
    head = NULL;
  }
  catch(...)
  {
    cerr << "Error in hal get chroms";
    head = NULL;
  }
  HAL_UNLOCK
  return head;
}

extern "C" char *halGetDna(int halHandle,
                           char* speciesName, char* chromName, 
                           hal_int_t start, hal_int_t end)
{
  HAL_LOCK
  char* dna = NULL;
  try
  {
    AlignmentConstPtr alignment = getExistingAlignment(halHandle, 0, true);
    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
      stringstream ss;
      ss << "Species with name " << speciesName << " not found in alignment "
         << "with handle " << halHandle;
      throw hal_exception(ss.str());
    }
    const Sequence* sequence = genome->getSequence(chromName);
    if (sequence == NULL)
    {
      stringstream ss;
      ss << "Chromosome with name " << chromName << " not found in species "
         << speciesName;
      throw hal_exception(ss.str());
    }    
    if (start > end || end > (hal_index_t)sequence->getSequenceLength())
    {
      stringstream ss;
      ss << "Specified range [" << start << "," << end << ") is invalid "
         << "for chromsome " << chromName << " in species " << speciesName
         << " which is of length " << sequence->getSequenceLength();
      throw hal_exception(ss.str());
    }
    
    string buffer;
    sequence->getSubString(buffer, start, end - start);
    dna = copyCString(buffer);
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    dna = NULL;
  }
  catch(...)
  {
    cerr << "Error in hal get dna";
    dna = NULL;
  }
  HAL_UNLOCK
  return dna;
}

void checkHandle(int handle)
{
  HandleMap::iterator mapIt = handleMap.find(handle);
  if (mapIt == handleMap.end())
  {
    stringstream ss;
    ss << "Handle " << handle << "not found in alignment map";
    throw hal_exception(ss.str());
  }
  if (mapIt->second.second.get() == NULL)
  {
    stringstream ss;
    ss << "Handle " << handle << "points to NULL alignment";
    throw hal_exception(ss.str());
  }
}

void checkGenomes(int halHandle, 
                  AlignmentConstPtr alignment, const string& qSpecies,
                  const string& tSpecies, const string& tChrom)
{
  const Genome* qGenome = alignment->openGenome(qSpecies);
  if (qGenome == NULL)
  {
    stringstream ss;
    ss << "Query species " << qSpecies << " not found in alignment with "
       << "handle " << halHandle;
    throw hal_exception(ss.str());
  }
  const Genome* tGenome = alignment->openGenome(tSpecies);
  if (tGenome == NULL)
  {
    stringstream ss;
    ss << "Reference species " << tSpecies << " not found in alignment with "
       << "handle " << halHandle;
    throw hal_exception(ss.str());
  }

  const Sequence* tSequence = tGenome->getSequence(tChrom);
  // cactus pipeline presently adds species name as prefix of 
  // sequence name.  check if this caused confusion
  string sequenceName;
  if (tSequence == NULL)
  {
    sequenceName = tGenome->getName();
    sequenceName += '.';
    sequenceName += tChrom;
    tSequence = tGenome->getSequence(sequenceName);
  }
  if (tSequence == NULL)
  {
    stringstream ss;
    ss << "Unable to locate sequence " << tChrom << "(or " << sequenceName 
       << ") in genome " << tSpecies;
    throw hal_exception(ss.str());
  }
}


AlignmentConstPtr getExistingAlignment(int handle, hal_size_t queryLength,
                                       bool needDNASequence)
{
  checkHandle(handle);
  HandleMap::iterator mapIt = handleMap.find(handle);
  return mapIt->second.second->getAlignment(queryLength, needDNASequence);
}

char* copyCString(const string& inString)
{
  char* outString = (char*)malloc(inString.length() + 1);
  strcpy(outString, inString.c_str());
  return outString;
}

hal_block_results_t* readBlocks(const Sequence* tSequence,
                                hal_index_t absStart, hal_index_t absEnd,
                                const Genome* qGenome, bool getSequenceString,
                                bool doDupes)
{
  const Genome* tGenome = tSequence->getGenome();
  string qGenomeName = qGenome->getName();
  hal_block_t* prev = NULL;
  BlockMapper blockMapper;
  blockMapper.init(tGenome, qGenome, absStart, absEnd, doDupes, 0, true);
  blockMapper.map();
  BlockMapper::MSSet paraSet;
  if (doDupes == true)
  {
    blockMapper.extractReferenceParalogies(paraSet);
  }
  BlockMapper::MSSet& segMap = blockMapper.getMap();
  vector<MappedSegmentConstPtr> fragments;

  hal_block_results_t* results = 
     (hal_block_results_t*)calloc(1, sizeof(hal_block_results_t));

  for (BlockMapper::MSSet::iterator segMapIt = segMap.begin();
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
    blockMapper.extractSegment(segMapIt, paraSet, fragments, &segMap);
    readBlock(cur, fragments, getSequenceString, qGenomeName);
    prev = cur;
  }

  results->targetDupeBlocks = processTargetDupes(blockMapper, paraSet);
  return results;
}

void readBlock(hal_block_t* cur,  
               vector<MappedSegmentConstPtr>& fragments, 
               bool getSequenceString, const string& genomeName)
{
  MappedSegmentConstPtr firstQuerySeg = fragments.front();
  MappedSegmentConstPtr lastQuerySeg = fragments.back();
  SlicedSegmentConstPtr firstRefSeg = firstQuerySeg->getSource();
  SlicedSegmentConstPtr lastRefSeg = lastQuerySeg->getSource();
  const Sequence* qSequence = firstQuerySeg->getSequence();
  const Sequence* tSequence = firstRefSeg->getSequence(); 
  assert(qSequence == lastQuerySeg->getSequence());
  assert(tSequence == lastRefSeg->getSequence());
  assert(firstRefSeg->getReversed() == false);
  assert(lastRefSeg->getReversed() == false);
  
  cur->next = NULL;

  string seqBuffer = qSequence->getName();
  string dnaBuffer;
  size_t prefix = 
     seqBuffer.find(genomeName + '.') != 0 ? 0 : genomeName.length() + 1;
  cur->qChrom = (char*)malloc(seqBuffer.length() + 1 - prefix);
  strcpy(cur->qChrom, seqBuffer.c_str() + prefix);

  cur->tStart = firstRefSeg->getStartPosition() - tSequence->getStartPosition();
  if (firstQuerySeg->getReversed() == false)
  {
    cur->qStart = firstQuerySeg->getStartPosition() - 
       qSequence->getStartPosition();
  }
  else
  {
    cur->qStart = lastQuerySeg->getEndPosition() - qSequence->getStartPosition();
  }
  assert(cur->tStart >= 0);
  assert(cur->qStart >= 0);

  assert(firstRefSeg->getLength() == firstQuerySeg->getLength());
  cur->size = lastRefSeg->getEndPosition() - firstRefSeg->getStartPosition() + 1;
  cur->strand = firstQuerySeg->getReversed() ? '-' : '+';
  cur->sequence = NULL;
  if (getSequenceString != 0)
  {
    qSequence->getSubString(dnaBuffer, cur->qStart, cur->size);
    if (cur->strand == '-')
    {
      reverseComplement(dnaBuffer);
    }
    cur->sequence = (char*)malloc(dnaBuffer.length() * sizeof(char));
    strcpy(cur->sequence, dnaBuffer.c_str());
  }
}

struct DupeIdLess { bool operator()(const hal_target_dupe_list_t* d1,
                                    const hal_target_dupe_list_t* d2) const {
  return d1->id < d2->id;
}};

hal_target_dupe_list_t* processTargetDupes(BlockMapper& blockMapper,
                                           BlockMapper::MSSet& paraSet)
{
  vector<hal_target_dupe_list_t*> tempList;
  vector<MappedSegmentConstPtr> fragments;
  BlockMapper::MSSet emptySet;

  // make a dupe list for each merged segment
  for (BlockMapper::MSSet::iterator segMapIt = paraSet.begin();
       segMapIt != paraSet.end(); ++segMapIt)
  {
    assert((*segMapIt)->getSource()->getReversed() == false);
    hal_target_dupe_list_t* cur = (hal_target_dupe_list_t*)calloc(
      1, sizeof(hal_target_dupe_list_t));
    blockMapper.extractSegment(segMapIt, emptySet, fragments, &paraSet);
    readTargetRange(cur, fragments);
    tempList.push_back(cur);
  }
  
  // sort based on query coordinate
  std::sort(tempList.begin(), tempList.end(), DupeIdLess());
  
  hal_target_dupe_list_t* head = NULL;
  hal_target_dupe_list_t* prev = NULL;
  int id = 0;
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
    while (j != tempList.end() && 
           ((*i)->id <= (*j)->id && (*i)->id + (*i)->tRange->size >= (*j)->id))
    {
      assert((*i)->next == NULL);
      assert((*j)->next == NULL);
      assert(j != i);
      k = j;
      ++k;
      (*j)->tRange->next = (*i)->tRange;
      (*i)->tRange = (*j)->tRange;
      (*j)->tRange = NULL;
      halFreeTargetDupeLists(*j);
      *j = NULL;
      j = k;
    }

    (*i)->id = id++;
    if (head == NULL)
    {
      head = *i;
    }
    else
    {
      prev->next = *i;
    }
    prev = *i;
    prev->next = NULL;
  }
  return head;
}

void readTargetRange(hal_target_dupe_list_t* cur,
                     vector<MappedSegmentConstPtr>& fragments)
{
  MappedSegmentConstPtr firstQuerySeg = fragments.front();
  MappedSegmentConstPtr lastQuerySeg = fragments.back();
  SlicedSegmentConstPtr firstRefSeg = firstQuerySeg->getSource();
  SlicedSegmentConstPtr lastRefSeg = lastQuerySeg->getSource();
  const Sequence* qSequence = firstQuerySeg->getSequence();
  const Sequence* tSequence = firstRefSeg->getSequence(); 
  assert(qSequence == lastQuerySeg->getSequence());
  assert(tSequence == lastRefSeg->getSequence());
  assert(firstRefSeg->getReversed() == false);
  assert(lastRefSeg->getReversed() == false);

  cur->tRange = (hal_target_range_t*)calloc(1, sizeof(hal_target_range_t));
   
  cur->tRange->tStart = firstRefSeg->getStartPosition() - 
     tSequence->getStartPosition();
  cur->tRange->size = lastRefSeg->getEndPosition() - 
     firstRefSeg->getStartPosition() + 1;

  // use query start as proxy for unique id right now
  if (firstQuerySeg->getReversed() == false)
  {
    cur->id = firstQuerySeg->getStartPosition() - 
       qSequence->getStartPosition();
  }
  else
  {
    cur->id = lastQuerySeg->getEndPosition() - qSequence->getStartPosition();
  }
}
