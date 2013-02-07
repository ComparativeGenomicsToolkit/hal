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

// We bend over backwards to cope with stone-age interface
typedef pair<AlignmentConstPtr, const Genome*> OpenGenome;
struct OGLess { 
   bool operator()(const OpenGenome& og1, const OpenGenome& og2) {
     return og1.first.get() < og2.first.get() || (
       og1.first.get() == og2.first.get() && og1.second < og2.second);
   }
};
  
static map<string, AlignmentConstPtr> halPathMap;
static map<int, OpenGenome> halHandleMap;
static map<OpenGenome, int, OGLess> halHandleReverseMap;

static block* readBlocks(BottomSegmentIteratorConstPtr bottom,
                         hal_index_t childIndex, hal_index_t absEnd,
                         int getSequenceString);
static void readBlock(block* cur, BottomSegmentIteratorConstPtr bottom,
                      TopSegmentIteratorConstPtr top, int getSequenceString,
                      const string& genomeName, string& seqBuffer, 
                      string& dnaBuffer);

extern "C" int halOpen(char *halFileName, char* qSpeciesName)
{
  int handle = -1;
  try
  {
    // if path has been opened before, find the alignment
    AlignmentConstPtr alignment;
    map<string, AlignmentConstPtr>::iterator pmi = 
       halPathMap.find(halFileName);
    if (pmi != halPathMap.end())
    {
      alignment = pmi->second;
    }
    else
    {
      alignment = hdf5AlignmentInstanceReadOnly();
      alignment->open(halFileName);
      halPathMap.insert(pair<string, AlignmentConstPtr>(
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
    map<OpenGenome, int, OGLess>::iterator hri = halHandleReverseMap.find(og);
    if (hri != halHandleReverseMap.end())
    {
      handle = hri->second;
    }
    else
    {
      // create new handle and return it.  updating both maps
      handle = 0;      
      if (halHandleMap.empty() == false)
      {
        handle = halHandleMap.rbegin()->first + 1;
      }
      halHandleMap.insert(pair<int, OpenGenome>(handle, og));
      halHandleReverseMap.insert(pair<OpenGenome, int>(og, handle));
    }
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    handle = -1;
  }
  
  return handle;
}

extern "C" int halClose(int handle)
{
  try
  {
    map<int, OpenGenome>::iterator hmi = halHandleMap.find(handle);
    if (hmi == halHandleMap.end())
    {
      stringstream ss;
      ss << "error closing handle " << handle << ": not found"; 
      throw hal_exception(ss.str());
    }
    OpenGenome og(hmi->second.first, hmi->second.second);
    AlignmentConstPtr alignment = og.first;
    // remove from handle maps
    halHandleMap.erase(hmi);
    halHandleReverseMap.erase(halHandleReverseMap.find(og));
    // close the genome
    alignment->closeGenome(og.second);

    bool last = true;
    for (hmi = halHandleMap.begin(); last && hmi != halHandleMap.end(); ++hmi)
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
      map<string, AlignmentConstPtr>::iterator hpi = halPathMap.begin();
      for (; hpi != halPathMap.end(); ++hpi)
      {
        AlignmentConstPtr alInMap = hpi->second;
        if (alInMap.get() == alignment.get())
        {
          alignment->close();
          halPathMap.erase(hpi);
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
                                                   int getSequenceString)
{
  try
  {
    map<int, OpenGenome>::iterator hmi = halHandleMap.find(halHandle);
    if (hmi == halHandleMap.end())
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
      sequenceName = tChrom;
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

    hal_index_t absStart = sequence->getStartPosition() + tStart;
    hal_index_t absEnd = sequence->getStartPosition() + tEnd;
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
    if (tEnd <= bottom->getEndPosition())
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
    string genomeName = top->getGenome()->getName();
    string seqBuffer, dnaBuffer;

    while (bottom->getStartPosition() < absEnd)
    {
      if (bottom->hasChild(childIndex))
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
        top->toChild(bottom, childIndex);
        readBlock(cur, bottom, top, getSequenceString, genomeName,
                  seqBuffer, dnaBuffer);        
        prev = cur;
        bottom->toRight(absEnd);
      }
    }    
    
    return head;
}

void readBlock(block* cur, BottomSegmentIteratorConstPtr bottom,
                      TopSegmentIteratorConstPtr top, int getSequenceString,
                      const string& genomeName, string& seqBuffer, 
                      string& dnaBuffer)
{
  const Sequence* tSequence = top->getSequence();
  cur->next = NULL;
  seqBuffer = tSequence->getName();
  size_t prefix = 
     seqBuffer.find(genomeName + '.') != 0 ? 0 : genomeName.length() + 1;
  cur->qChrom = (char*)malloc(seqBuffer.length() + 1 - prefix);
  strcpy(cur->qChrom, seqBuffer.c_str() + prefix);
  cur->tStart = bottom->getStartPosition();
  cur->qStart = top->getStartPosition();
  if (top->getReversed() == true)
  {
    cur->qStart = tSequence->getSequenceLength() - cur->qStart;
  }
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
