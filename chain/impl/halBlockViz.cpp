/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <map>
#include <sstream>
#include "hal.h"
#include "halChain.h"
#include "halBlockViz.h"
#include "halBlockMapper.h"

using namespace std;
using namespace hal;

typedef map<int, AlignmentConstPtr> HandleMap;
static int MAX_THREAD_ID = -1;
static HandleMap* HANDLE_MAP_ARRAY = NULL;

static void checkThreadID(int threadID);
static void checkHandle(int handle, int threadID);
static AlignmentConstPtr getExistingAlignment(int handle, int threadID);
static char* copyCString(const string& inString);

static hal_block_t* readBlocks(const Sequence* tSequence,
                               hal_index_t absStart, hal_index_t absEnd,
                               const Genome* qGenome, bool getSequenceString,
                               bool doDupes);

static void readBlock(hal_block_t* cur, SegmentIteratorConstPtr refSeg,
                      SegmentIteratorConstPtr querySeg, bool getSequenceString,
                      const string& genomeName);



extern "C" int halInit(int maxThreadID)
{
  try
  {
    halExit();
    MAX_THREAD_ID = maxThreadID;
    HANDLE_MAP_ARRAY = new HandleMap[MAX_THREAD_ID + 1];
  }
  catch(...)
  {
    return -1;
  }
  return 0;
}

extern "C" int halExit()
{
  try
  {
    for (int i = 0; i <= MAX_THREAD_ID; ++i)
    {
      for (HandleMap::iterator j = HANDLE_MAP_ARRAY[i].begin(); 
           j != HANDLE_MAP_ARRAY[i].end(); ++j)
      {
        j->second->close();
      }
    }
    delete [] HANDLE_MAP_ARRAY;
    HANDLE_MAP_ARRAY = NULL;
    MAX_THREAD_ID = -1;
  }
  catch(...)
  {
    return -1;
  }
  return 0;
}

extern "C" int halOpen(char* halFileName, int threadID)
{
  int handle = -1;
  try
  {
    checkThreadID(threadID);
    HandleMap& handleMap = HANDLE_MAP_ARRAY[threadID];
    HandleMap::reverse_iterator mapIt = handleMap.rbegin();
    if (mapIt == handleMap.rend())
    {
      handle = 0;
    }
    else
    {
      handle = mapIt->first + 1;
    }
    
    AlignmentConstPtr alignment = hdf5AlignmentInstanceReadOnly();
    alignment->open(halFileName);
    handleMap.insert(pair<int, AlignmentConstPtr>(handle, alignment));
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    handle = -1;
  }

  return handle;
}

extern "C" int halClose(int handle, int threadID)
{
  try
  {
    checkHandle(handle, threadID);
    HandleMap& handleMap = HANDLE_MAP_ARRAY[threadID];
    HandleMap::iterator mapIt = handleMap.find(handle);
    if (mapIt == handleMap.end())
    {
      stringstream ss;
      ss << "error closing handle " << handle << ": not found";
      throw hal_exception(ss.str());
    }
    mapIt->second->close();
    handleMap.erase(mapIt);
  }
   catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return -1;
  }

  return 0;
}

extern "C" void halFreeBlocks(hal_block_t* head)
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

extern "C" struct hal_block_t *halGetBlocksInTargetRange(int halHandle,
                                                         int threadID,
                                                         char* qSpecies,
                                                         char* tSpecies,
                                                         char* tChrom,
                                                         int tStart, int tEnd,
                                                         int getSequenceString,
                                                         int doDupes)
{
  try
  {
    AlignmentConstPtr alignment = getExistingAlignment(halHandle, threadID);
    const Genome* qGenome = alignment->openGenome(qSpecies);
    if (qGenome == NULL)
    {
      stringstream ss;
      ss << "Query species " << qSpecies << " not found in alignment with "
         << "handle " << halHandle << " in thread " << threadID;
      throw hal_exception(ss.str());
    }
    const Genome* tGenome = alignment->openGenome(tSpecies);
    if (tGenome == NULL)
    {
      stringstream ss;
      ss << "Reference species " << tSpecies << " not found in alignment with "
         << "handle " << halHandle << " in thread " << threadID;
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
    if (tSequence == NULL || 
        (hal_size_t)tStart >= tSequence->getSequenceLength() ||
        (hal_size_t)tEnd > tSequence->getSequenceLength())
    {
      stringstream ss;
      ss << "Unable to locate sequence " << tChrom << "(or " << sequenceName 
         << ") in genome " << tSpecies;
      throw hal_exception(ss.str());
    }

    hal_index_t myEnd = tEnd > 0 ? tEnd : tSequence->getSequenceLength();
    hal_index_t absStart = tSequence->getStartPosition() + tStart;
    hal_index_t absEnd = tSequence->getStartPosition() + myEnd - 1;
    if (absStart > absEnd)
    {
      throw hal_exception("Invalid range");
    }

    hal_block_t* head = readBlocks(tSequence, absStart, absEnd, qGenome,
                                   getSequenceString != 0, doDupes != 0);

    return head;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return NULL;
}

extern "C" struct hal_species_t *halGetSpecies(int halHandle, int threadID)
{
  try
  {
    AlignmentConstPtr alignment = getExistingAlignment(halHandle, threadID);
    hal_species_t* head = NULL;
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
    return head;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return NULL;
}

extern "C" struct hal_chromosome_t *halGetChroms(int halHandle, int threadID,
                                                 char* speciesName)
{
  try
  {
    AlignmentConstPtr alignment = getExistingAlignment(halHandle, threadID);
    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
      stringstream ss;
      ss << "Species with name " << speciesName << " not found in alignment "
         << "with handle " << halHandle << " and threadID " << threadID;
      throw hal_exception(ss.str());
    }

    hal_chromosome_t* head = NULL;
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
    return head;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return NULL;
}

extern "C" char *halGetDna(int halHandle, int threadID,
                           char* speciesName, char* chromName, 
                           int start, int end)
{
  try
  {
    AlignmentConstPtr alignment = getExistingAlignment(halHandle, threadID);
    const Genome* genome = alignment->openGenome(speciesName);
    if (genome == NULL)
    {
      stringstream ss;
      ss << "Species with name " << speciesName << " not found in alignment "
         << "with handle " << halHandle << " and threadID " << threadID;
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
    return copyCString(buffer);
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return NULL;
}

void checkThreadID(int threadID)
{
  if (HANDLE_MAP_ARRAY == NULL)
  {
    stringstream ss;
    ss << "halInit() must be called before using any other methods in this API";
    throw hal_exception(ss.str());
  }
  if (threadID > MAX_THREAD_ID)
  {
    stringstream ss;
    ss << "Max thread ID set by halInit() is " <<MAX_THREAD_ID
       << ":  Input threadID " << threadID << " is out of range.";
    throw hal_exception(ss.str());
  }
}

void checkHandle(int handle, int threadID)
{
  checkThreadID(threadID);
  HandleMap& handleMap = HANDLE_MAP_ARRAY[threadID];
  HandleMap::iterator mapIt = handleMap.find(handle);
  if (mapIt == handleMap.end())
  {
    stringstream ss;
    ss << "Handle " << handle << "not found in alignment map for thread id "
       << threadID;
    throw hal_exception(ss.str());
  }
  if (mapIt->second.get() == NULL)
  {
    stringstream ss;
    ss << "Handle " << handle << "points to NULL alignment for thread id "
       << threadID;
    throw hal_exception(ss.str());
  }
}

AlignmentConstPtr getExistingAlignment(int handle, int threadID)
{
  checkHandle(handle, threadID);
  HandleMap& handleMap = HANDLE_MAP_ARRAY[threadID];
  HandleMap::iterator mapIt = handleMap.find(handle);
  return mapIt->second;
}

char* copyCString(const string& inString)
{
  char* outString = (char*)malloc(inString.length() + 1);
  strcpy(outString, inString.c_str());
  return outString;
}

hal_block_t* readBlocks(const Sequence* tSequence,
                        hal_index_t absStart, hal_index_t absEnd,
                        const Genome* qGenome, bool getSequenceString,
                        bool doDupes)
{
  const Genome* tGenome = tSequence->getGenome();
  string qGenomeName = qGenome->getName();
  hal_block_t* head = NULL;
  hal_block_t* prev = NULL;
  BlockMapper blockMapper;
  blockMapper.init(tGenome, qGenome, absStart, absEnd, false);
  blockMapper.map();
  const BlockMapper::SegMap& segMap = blockMapper.getMap();
  for (BlockMapper::SegMap::const_iterator segMapIt = segMap.begin();
       segMapIt != segMap.end(); ++segMapIt)
  {
    SegmentIteratorConstPtr refSeg = segMapIt->first;
    BlockMapper::SegSet* segSet = segMapIt->second;
    assert(doDupes || segSet->size() == 1);
    for (BlockMapper::SegSet::const_iterator segIt = segSet->begin();
         segIt != segSet->end(); ++segIt)
    {
      hal_block_t* cur = (hal_block_t*)calloc(1, sizeof(hal_block_t));
      if (head == NULL)
      {
        head = cur;
      }
      else
      {
        prev->next = cur;
      }
      readBlock(cur, refSeg, *segIt, getSequenceString, qGenomeName);
      prev = cur;
    }    
  }
  return head;
}

void readBlock(hal_block_t* cur, SegmentIteratorConstPtr refSeg,
               SegmentIteratorConstPtr querySeg, bool getSequenceString,
               const string& genomeName)
{
  const Sequence* qSequence = querySeg->getSequence();
  const Sequence* tSequence = refSeg->getSequence(); 
  cur->next = NULL;

  string seqBuffer = qSequence->getName();
  string dnaBuffer;
  size_t prefix = 
     seqBuffer.find(genomeName + '.') != 0 ? 0 : genomeName.length() + 1;
  cur->qChrom = (char*)malloc(seqBuffer.length() + 1 - prefix);
  strcpy(cur->qChrom, seqBuffer.c_str() + prefix);

  cur->tStart = refSeg->getStartPosition() - tSequence->getStartPosition();
  cur->qStart = querySeg->getStartPosition() - qSequence->getStartPosition();
  assert(cur->tStart >= 0);
  assert(cur->qStart >= 0);

  assert(refSeg->getLength() == querySeg->getLength());
  cur->size = refSeg->getLength();
  cur->strand = querySeg->getReversed() ? '-' : '+';
  cur->sequence = NULL;
  if (getSequenceString != 0)
  {
    querySeg->getString(dnaBuffer);
    cur->sequence = (char*)malloc(dnaBuffer.length() + 1);
    strcpy(cur->sequence, dnaBuffer.c_str());
  }
}

