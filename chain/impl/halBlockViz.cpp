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

using namespace std;
using namespace hal;

static const int MAX_CLIENT_THREADS = 10000;

// We bend over backwards to cope with stone-age interface
typedef pair<AlignmentConstPtr, const Genome*> OpenGenome;
struct OGLess { 
   bool operator()(const OpenGenome& og1, const OpenGenome& og2) {
     return og1.first.get() < og2.first.get() || (
       og1.first.get() == og2.first.get() && og1.second < og2.second);
   }
};
  
static map<string, AlignmentConstPtr> halPathMap[MAX_CLIENT_THREADS];
static map<int, OpenGenome> halHandleMap[MAX_CLIENT_THREADS];
static map<OpenGenome, int, OGLess> halHandleReverseMap[MAX_CLIENT_THREADS];

static block* readBlocks(BottomSegmentIteratorConstPtr bottom,
                         hal_index_t childIndex, hal_index_t absEnd,
                         int getSequenceString);
static void readBlock(block* cur, BottomSegmentIteratorConstPtr bottom,
                      TopSegmentIteratorConstPtr top, int getSequenceString,
                      const string& genomeName, string& seqBuffer, 
                      string& dnaBuffer);
static bool toAdjacent(BottomSegmentIteratorConstPtr bottom, 
                       hal_index_t childIndex,
                       TopSegmentIteratorConstPtr top,
                       bool left);


extern "C" int halOpen(char *halFileName, char* qSpeciesName, int threadID)
{
  int handle = -1;
  try
  {
    if (threadID >= MAX_CLIENT_THREADS)
    {
      stringstream ss;
      ss << "Thread ID " << threadID << " exceeds the hardcoded limit "
         << MAX_CLIENT_THREADS << " Please number threadIDs from 0...N or "
         << "increase MAX_CLIENT_THREADS in halBlockViz.cpp";
      throw hal_exception(ss.str());
    }
    
    // if path has been opened before, find the alignment
    AlignmentConstPtr alignment;
    map<string, AlignmentConstPtr>::iterator pmi =
       halPathMap[threadID].find(halFileName);
    if (pmi != halPathMap[threadID].end())
    {
      alignment = pmi->second;
    }
    else
    {
      alignment = hdf5AlignmentInstanceReadOnly();
      alignment->open(halFileName);
      halPathMap[threadID].insert(pair<string, AlignmentConstPtr>(
                                    qSpeciesName, alignment));
    }

    const Genome* genome = alignment->openGenome(qSpeciesName);
    if (genome == NULL)
    {
      throw hal_exception(string(qSpeciesName) + " not found");
    }
    OpenGenome og(alignment, genome);
    // if the path and genome have been open before, return the existing
    // handle
    map<OpenGenome, int, OGLess>::iterator hri = 
       halHandleReverseMap[threadID].find(og);
    if (hri != halHandleReverseMap[threadID].end())
    {
      handle = hri->second;
    }
    else
    {
      // create new handle and return it.  updating both maps
      handle = 0;
      if (halHandleMap[threadID].empty() == false)
      {
        handle = halHandleMap[threadID].rbegin()->first + 1;
      }
      halHandleMap[threadID].insert(pair<int, OpenGenome>(handle, og));
      halHandleReverseMap[threadID].insert(pair<OpenGenome, int>(og, handle));
    }
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
    if (threadID >= MAX_CLIENT_THREADS)
    {
      stringstream ss;
      ss << "Thread ID " << threadID << " exceeds the hardcoded limit "
         << MAX_CLIENT_THREADS << " Please number threadIDs from 0...N or "
         << "increase MAX_CLIENT_THREADS in halBlockViz.cpp";
      throw hal_exception(ss.str());
    }

    map<int, OpenGenome>::iterator hmi = halHandleMap[threadID].find(handle);
    if (hmi == halHandleMap[threadID].end())
    {
      stringstream ss;
      ss << "error closing handle " << handle << ": not found";
      throw hal_exception(ss.str());
    }
    OpenGenome og(hmi->second.first, hmi->second.second);
    AlignmentConstPtr alignment = og.first;
    // remove from handle maps
    halHandleMap[threadID].erase(hmi);
    halHandleReverseMap[threadID].erase(halHandleReverseMap[threadID].find(og));
    // close the genome
    alignment->closeGenome(og.second);

    bool last = true;
    for (hmi = halHandleMap[threadID].begin(); 
         last && hmi != halHandleMap[threadID].end(); ++hmi)
    {
      AlignmentConstPtr alInMap = hmi->second.first;
      if (alInMap.get() == alignment.get())
      {
        last = false;
      }
    }

    // this was the last reference to the alignment. so we close it
    if (last == true)
    {
      map<string, AlignmentConstPtr>::iterator hpi = 
         halPathMap[threadID].begin();
      for (; hpi != halPathMap[threadID].end(); ++hpi)
      {
        AlignmentConstPtr alInMap = hpi->second;
        if (alInMap.get() == alignment.get())
        {
          alignment->close();
          halPathMap[threadID].erase(hpi);
          break;
        }
      }
    }
  }
   catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return -1;
  }


  return 0;
}

extern "C" void halFreeBlocks(block* head)
{
  while (head != NULL)
  {
    block* next = head->next;
    free(head->sequence);
    free(head->qChrom);
    free(head);
    head = next;
  }
}

extern "C" struct block *halGetBlocksInTargetRange(int halHandle,
                                                   char* tChrom,
                                                   int tStart, int tEnd,
                                                   int getSequenceString,
                                                   int threadID)
{
  try
  {
    map<int, OpenGenome>::iterator hmi = halHandleMap[threadID].find(halHandle);
    if (hmi == halHandleMap[threadID].end())
    {
      stringstream ss;
      ss << "Handle " << halHandle << " not found";
      throw hal_exception(ss.str());
    }
    AlignmentConstPtr alignment = hmi->second.first;
    const Genome* genome = hmi->second.second;
    const Genome* parent = genome->getParent();
    if (parent == NULL)
    {
      stringstream ss;
      ss << "Open genome " << genome->getName() << " is root in HAL tree"
         " and therefore has no reference and cannot be used to get blocks";
      throw hal_exception(ss.str());
    }

    const Sequence* sequence = parent->getSequence(tChrom);
    // cactus pipeline presently adds species name as prefix of 
    // sequence name.  check if this caused confusion
    string sequenceName;
    if (sequence == NULL)
    {
      sequenceName = parent->getName();
      sequenceName += '.';
      sequenceName += tChrom;
      sequence = parent->getSequence(sequenceName);
    }
    if (sequence == NULL || 
        (hal_size_t)tStart >= sequence->getSequenceLength() ||
        (hal_size_t)tEnd > sequence->getSequenceLength())
    {
      stringstream ss;
      ss << "Unable to locate sequence " << tChrom << "(or " << sequenceName 
         << ") in genome " << parent->getName();
      throw hal_exception(ss.str());
    }

    hal_index_t myEnd = tEnd > 0 ? tEnd : sequence->getSequenceLength();
    hal_index_t absStart = sequence->getStartPosition() + tStart;
    hal_index_t absEnd = sequence->getStartPosition() + myEnd;
    hal_index_t childIndex = parent->getChildIndex(genome);
    if (absStart >= absEnd)
    {
      throw hal_exception("Invalid range");
    }
    
    BottomSegmentIteratorConstPtr bottom = 
       parent->getBottomSegmentIterator(childIndex);
    bottom->toSite(absStart, false);
    hal_offset_t startOffset = absStart - bottom->getStartPosition();
    hal_offset_t endOffset = 0;
    if (absEnd <= bottom->getEndPosition())
    {
      endOffset = bottom->getEndPosition() - absEnd + 1;
    }

    bottom->slice(startOffset, endOffset);
    assert(bottom->getStartPosition() == absStart);
    assert(bottom->getEndPosition() <= absEnd);

    if (bottom->getLength() == 0)
    {
      throw  hal_exception("Error generating blocks.  "
                           "Invalid range specified?");
    }

    block* head = readBlocks(bottom, childIndex, absEnd, getSequenceString);
    return head;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return NULL;
}

block* readBlocks(BottomSegmentIteratorConstPtr bottom, hal_index_t childIndex,
                  hal_index_t absEnd, int getSequenceString)
{
  block* head = NULL;
  block* prev = NULL;

  TopSegmentIteratorConstPtr top = 
     bottom->getGenome()->getChild(childIndex)->getTopSegmentIterator();
  BottomSegmentIteratorConstPtr bottom2 = 
     bottom->getGenome()->getBottomSegmentIterator();
  string genomeName = top->getGenome()->getName();
  string seqBuffer, dnaBuffer;
  hal_index_t lastIndex = 
     (hal_index_t)bottom->getGenome()->getNumBottomSegments();
  std::set<hal_index_t> addedSet;
  while (bottom->getArrayIndex() < lastIndex &&
         bottom->getStartPosition() < absEnd - 1)
  {
    if (bottom->hasChild(childIndex))
    {
      top->toChild(bottom, childIndex);
      if (addedSet.find(top->getStartPosition()) == addedSet.end())
      {
        block* cur = (block*)malloc(sizeof(block));
        if (head == NULL)
        {
          head = cur;
        }
        else
        {
          prev->next = cur;
        }        
        readBlock(cur, bottom, top, getSequenceString, genomeName,
                  seqBuffer, dnaBuffer);
        addedSet.insert(top->getStartPosition());
       
        // insert next aligned block to left of current block
        // in array
        bool res = toAdjacent(bottom, childIndex, top, true);
        if (res && addedSet.find(top->getStartPosition()) == addedSet.end())
        {
          bottom2->toParent(top);
          block* left = (block*)malloc(sizeof(block));
          readBlock(left, bottom2, top, getSequenceString, genomeName,
                    seqBuffer, dnaBuffer);
          addedSet.insert(top->getStartPosition());
          if (prev == NULL)
          {
            head = left;
          }
          else
          {
            prev->next = left;
          }
          left->next = cur;
        }
        
        // insert next aligned block to right of current block
        // in array
        res = toAdjacent(bottom, childIndex, top, false);
        if (res && addedSet.find(top->getStartPosition()) == addedSet.end())
        {
          bottom2->toParent(top);
          block* right = (block*)malloc(sizeof(block));
          readBlock(right, bottom2, top, getSequenceString, genomeName,
                    seqBuffer, dnaBuffer);
          addedSet.insert(top->getStartPosition());
          cur->next = right;
          cur = right;
        }
        
        prev = cur;
      }
    }
    bottom->toRight(absEnd - 1);
  }    
    
  return head;
}

void readBlock(block* cur, BottomSegmentIteratorConstPtr bottom,
               TopSegmentIteratorConstPtr top, int getSequenceString,
               const string& genomeName, string& seqBuffer, 
               string& dnaBuffer)
{
  const Sequence* qSequence = top->getSequence();
  const Sequence* tSequence = bottom->getSequence(); 
  cur->next = NULL;

  seqBuffer = qSequence->getName();
  size_t prefix = 
     seqBuffer.find(genomeName + '.') != 0 ? 0 : genomeName.length() + 1;
  cur->qChrom = (char*)malloc(seqBuffer.length() + 1 - prefix);
  strcpy(cur->qChrom, seqBuffer.c_str() + prefix);

  cur->tStart = bottom->getStartPosition() - tSequence->getStartPosition();
  cur->qStart = top->getStartPosition() - qSequence->getStartPosition();
  assert(cur->tStart >= 0);
  assert(cur->qStart >= 0);

  assert(bottom->getLength() == top->getLength());
  cur->size = bottom->getLength();
  cur->strand = top->getReversed() ? '-' : '+';
  cur->sequence = NULL;
  if (getSequenceString != 0)
  {
    top->getString(dnaBuffer);
    cur->sequence = (char*)malloc(dnaBuffer.length() + 1);
    strcpy(cur->sequence, dnaBuffer.c_str());
  }
}

// next segment in query sequence that aligns to reference.  
// return false if doesnt exist.
bool toAdjacent(BottomSegmentIteratorConstPtr bottom, 
                hal_index_t childIndex,
                TopSegmentIteratorConstPtr top,
                bool left)
{
  top->toChild(bottom, childIndex);
  const Sequence* seq = top->getSequence();
  hal_index_t N = (hal_index_t)top->getGenome()->getNumTopSegments();
  while (true)
  {
    if ((left && !top->getReversed()) ||
        (!left && top->getReversed()))
    {
      top->toLeft();
    }
    else
    {
      top->toRight();
    }
    if (top->getArrayIndex() < 0 || top->getArrayIndex() >= N ||
        top->getSequence() != seq)
    {
      break;
    }
    if (top->hasParent())
    {
      return true;
    }    
  }

  return false;
}
