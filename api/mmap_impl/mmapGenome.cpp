#include "mmapGenome.h"
#include "mmapTopSegment.h"
#include "mmapBottomSegment.h"
#include "mmapSequence.h"
#include "mmapSequenceIterator.h"
#include "mmapDNAIterator.h"
#include "defaultTopSegmentIterator.h"
#include "defaultBottomSegmentIterator.h"
#include "defaultColumnIterator.h"
#include "defaultRearrangement.h"
#include "defaultGappedTopSegmentIterator.h"
#include "defaultGappedBottomSegmentIterator.h"
using namespace hal;
using namespace std;

void MMapGenome::setDimensions(
  const vector<Sequence::Info>& sequenceDimensions,
  bool storeDNAArrays)
{
    hal_size_t totalSequenceLength = 0;
   hal_size_t maxName = 0;

    // Copy segment dimensions to use the external interface
    vector<Sequence::UpdateInfo> topDimensions;
    topDimensions.reserve(sequenceDimensions.size());
    vector<Sequence::UpdateInfo> bottomDimensions;
    bottomDimensions.reserve(sequenceDimensions.size());

    // Compute summary info from the list of sequence Dimensions
    for (vector<Sequence::Info>::const_iterator i = sequenceDimensions.begin();
         i != sequenceDimensions.end(); 
         ++i)
    {
        totalSequenceLength += i->_length;
        maxName = max(static_cast<hal_size_t>(i->_name.length()), maxName);
        topDimensions.push_back(
            Sequence::UpdateInfo(i->_name, i->_numTopSegments));
        bottomDimensions.push_back(
            Sequence::UpdateInfo(i->_name, i->_numBottomSegments));
    }

    // Write the new DNA/sequence information.
    _data->_totalSequenceLength = totalSequenceLength;
    _data->_dnaOffset = _alignment->allocateNewArray(sizeof(char) * totalSequenceLength);
    // FIXME: this needs to replaced, for two reasons: a) alignment
    // and b) Sequence class needs to carry around more information,
    // like the number of top/bottom segments for some reason.
    size_t sequenceElementSize = maxName + 1 + sizeof(hal_size_t);
    _data->_sequencesOffset = _alignment->allocateNewArray(sequenceElementSize * sequenceDimensions.size());
    hal_size_t startPos = 0;
    for (size_t i = 0;
         i < sequenceDimensions.size();
         ++i)
    {
        setSequenceElement(i, sequenceDimensions[i]._name, startPos);
        startPos += sequenceDimensions[i]._length;
    }

    // Write the new segment data.
    updateTopDimensions(topDimensions);
    updateBottomDimensions(bottomDimensions);
}

void MMapGenome::setSequenceElement(size_t i, const std::string &name, hal_size_t startPos) {
    // TODO
}

void MMapGenome::updateTopDimensions(
  const vector<Sequence::UpdateInfo>& topDimensions)
{
    // TODO
}

void MMapGenome::updateBottomDimensions(
  const vector<Sequence::UpdateInfo>& bottomDimensions)
{
    // TODO
}

hal_size_t MMapGenome::getNumSequences() const
{
  return _data->_numSequences;
}

Sequence* MMapGenome::getSequence(const string& name)
{
  loadSequenceNameCache();
  Sequence* sequence = NULL;
  map<string, MMapSequence*>::iterator mapIt = _sequenceNameCache.find(name);
  if (mapIt != _sequenceNameCache.end())
  {
    sequence = mapIt->second;
  }
  return sequence;
}

const Sequence* MMapGenome::getSequence(const string& name) const
{
  loadSequenceNameCache();
  const Sequence* sequence = NULL;
  map<string, MMapSequence*>::const_iterator mapIt = 
     _sequenceNameCache.find(name);
  if (mapIt != _sequenceNameCache.end())
  {
    sequence = mapIt->second;
  }
  return sequence;
}

Sequence* MMapGenome::getSequenceBySite(hal_size_t position)
{
  loadSequencePosCache();
  map<hal_size_t, MMapSequence*>::iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= (hal_size_t)i->second->getStartPosition())
    {
      assert(position < i->second->getStartPosition() +
             i->second->getSequenceLength());
      return i->second;
    }
  }
  return NULL;
}

const Sequence* MMapGenome::getSequenceBySite(hal_size_t position) const
{
  loadSequencePosCache();
  map<hal_size_t, MMapSequence*>::const_iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= (hal_size_t)i->second->getStartPosition())
    {
      assert(position < i->second->getStartPosition() +
             i->second->getSequenceLength());
      return i->second;
    }
  }
  return NULL;
}

SequenceIteratorPtr MMapGenome::getSequenceIterator(
  hal_index_t position)
{
  assert(position <= (hal_index_t)_sequenceNameArray.getSize());
  MMapSequenceIterator* newIt = new MMapSequenceIterator(this, position);
  return SequenceIteratorPtr(newIt);
}

SequenceIteratorConstPtr MMapGenome::getSequenceIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)_sequenceNameArray.getSize());
  // genome effectively gets re-consted when returned in the
  // const iterator.  just save doubling up code.
  MMapSequenceIterator* newIt = new MMapSequenceIterator(
    const_cast<MMapGenome*>(this), position);
  return SequenceIteratorConstPtr(newIt);
}

SequenceIteratorConstPtr MMapGenome::getSequenceEndIterator() const
{
  return getSequenceIterator(getNumSequences());
}

MetaData* MMapGenome::getMetaData()
{
    // TODO
    throw hal_exception("All metadata functions currently unimplemented");
}

const MetaData* MMapGenome::getMetaData() const
{
    // TODO
    throw hal_exception("All metadata functions currently unimplemented");
}

Genome* MMapGenome::getParent()
{
  return _alignment->openGenome(_alignment->getParentName(getName()));
}

const Genome* MMapGenome::getParent() const
{
  return _alignment->openGenome(_alignment->getParentName(getName()));
}

Genome* MMapGenome::getChild(hal_size_t childIdx)
{
  assert(childIdx < _data->_numChildren);
  vector<string> childNames = _alignment->getChildNames(getName());
  assert(childNames.size() > childIdx);
  return _alignment->openGenome(childNames.at(childIdx));
}

const Genome* MMapGenome::getChild(hal_size_t childIdx) const
{
  assert(childIdx < _data->_numChildren);
  vector<string> childNames = _alignment->getChildNames(getName());
  assert(childNames.size() > childIdx);
  return _alignment->openGenome(childNames.at(childIdx));
}

hal_size_t MMapGenome::getNumChildren() const
{
  return _data->_numChildren;
}

hal_index_t MMapGenome::getChildIndex(const Genome* child) const
{
  string childName = child->getName();
  vector<string> childNames = _alignment->getChildNames(getName());
  for (hal_size_t i = 0; i < childNames.size(); ++i)
  {
    if (childNames[i] == childName)
    {
      return i;
    }
  }
  return NULL_INDEX;
}

bool MMapGenome::containsDNAArray() const
{
  // FIXME: this will cause issues when there really isn't any DNA to
  // show, but I"m not sure we really want to support those use cases
  // anymore.
  return true;
}

const Alignment* MMapGenome::getAlignment() const
{
  return _alignment;
}

// SEGMENTED SEQUENCE INTERFACE

const string& MMapGenome::getName() const
{
    // TODO
    throw hal_exception("unimplemented");
}

hal_size_t MMapGenome::getSequenceLength() const
{
  return _data->_totalSequenceLength;
}

hal_size_t MMapGenome::getNumTopSegments() const
{
    return _data->_numTopSegments;
}

hal_size_t MMapGenome::getNumBottomSegments() const
{
  return _data->_numBottomSegments;
}

TopSegmentIteratorPtr MMapGenome::getTopSegmentIterator(hal_index_t segmentIndex)
{
  assert(position <= (hal_index_t)getNumTopSegments());
  MMapTopSegment *newSeg = new MMapTopSegment(this, segmentIndex);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultTopSegmentIterator *newIt = new DefaultTopSegmentIterator(newSeg);
  return TopSegmentIteratorPtr(newIt);
}

TopSegmentIteratorConstPtr MMapGenome::getTopSegmentIterator(hal_index_t segmentIndex) const
{
  assert(position <= (hal_index_t)getNumTopSegments());
  MMapGenome *genome = const_cast<MMapGenome*>(this);
  MMapTopSegment *newSeg = new MMapTopSegment(genome, segmentIndex);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultTopSegmentIterator* newIt = new DefaultTopSegmentIterator(newSeg);
  return TopSegmentIteratorConstPtr(newIt);
}

TopSegmentIteratorConstPtr MMapGenome::getTopSegmentEndIterator() const
{
  return getTopSegmentIterator(getNumTopSegments());
}

BottomSegmentIteratorPtr MMapGenome::getBottomSegmentIterator(hal_index_t segmentIndex)
{
  assert(position <= (hal_index_t)getNumBottomSegments());
  MMapBottomSegment* newSeg = new MMapBottomSegment(this, segmentIndex);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultBottomSegmentIterator* newIt = new DefaultBottomSegmentIterator(newSeg);
  return BottomSegmentIteratorPtr(newIt);
}

BottomSegmentIteratorConstPtr MMapGenome::getBottomSegmentIterator(hal_index_t segmentIndex) const
{
  assert(position <= (hal_index_t)getNumBottomSegments());
  MMapGenome* genome = const_cast<MMapGenome*>(this);
  MMapBottomSegment* newSeg = new MMapBottomSegment(genome, segmentIndex);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultBottomSegmentIterator* newIt = new DefaultBottomSegmentIterator(newSeg);
  return BottomSegmentIteratorConstPtr(newIt);
}

BottomSegmentIteratorConstPtr MMapGenome::getBottomSegmentEndIterator() const
{
  return getBottomSegmentIterator(getNumBottomSegments());
}
   
DNAIteratorPtr MMapGenome::getDNAIterator(hal_index_t position)
{
  assert(position / 2 <= (hal_index_t)_dnaArray.getSize());
  MMapDNAIterator* newIt = new MMapDNAIterator(this, position);
  return DNAIteratorPtr(newIt);
}

DNAIteratorConstPtr MMapGenome::getDNAIterator(hal_index_t position) const
{
  assert(_dnaArray.getSize() == 0 || 
         position / 2 <= (hal_index_t)_dnaArray.getSize());
  MMapGenome* genome = const_cast<MMapGenome*>(this);
  const MMapDNAIterator* newIt = new MMapDNAIterator(genome, position);
  return DNAIteratorConstPtr(newIt);
}

DNAIteratorConstPtr MMapGenome::getDNAEndIterator() const
{
  return getDNAIterator(getSequenceLength());
}

ColumnIteratorConstPtr MMapGenome::getColumnIterator(
  const set<const Genome*>* targets, hal_size_t maxInsertLength, 
  hal_index_t position, hal_index_t lastPosition, bool noDupes,
  bool noAncestors, bool reverseStrand, bool unique, bool onlyOrthologs) const
{
  hal_index_t lastIdx = lastPosition;
  if (lastPosition == NULL_INDEX)
  {
    lastIdx = (hal_index_t)(getSequenceLength() - 1);
  }
  if (position < 0 || 
      lastPosition >= (hal_index_t)(getSequenceLength()))
  {
    stringstream ss;
    ss << "MMapGenome::getColumnIterator: input indices "
       << "(" << position << ", " << lastPosition << ") out of bounds";
    throw hal_exception(ss.str());
  }
  const DefaultColumnIterator* newIt = 
     new DefaultColumnIterator(this, targets, position, lastIdx, 
                               maxInsertLength, noDupes, noAncestors,
                               reverseStrand, unique, onlyOrthologs);
  return ColumnIteratorConstPtr(newIt);
}

void MMapGenome::getString(string& outString) const
{
  getSubString(outString, 0, getSequenceLength());
}

void MMapGenome::setString(const string& inString)
{
  setSubString(inString, 0, getSequenceLength());
}

void MMapGenome::getSubString(string& outString, hal_size_t start,
                              hal_size_t length) const
{
  outString.resize(length);
  MMapDNAIterator dnaIt(const_cast<MMapGenome*>(this), start);
  dnaIt.readString(outString, length);
}

void MMapGenome::setSubString(const string& inString, 
                              hal_size_t start,
                              hal_size_t length)
{
  if (length != inString.length())
  {
    throw hal_exception(string("setString: input string has differnt") +
                               "length from target string in genome");
  }
  MMapDNAIterator dnaIt(this, start);
  dnaIt.writeString(inString, length);
}

RearrangementPtr MMapGenome::getRearrangement(hal_index_t position,
                                              hal_size_t gapLengthThreshold,
                                              double nThreshold,
                                              bool atomic) const
{
  assert(position >= 0 && position < (hal_index_t)getNumTopSegments());
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(position);  
  DefaultRearrangement* rea = new DefaultRearrangement(this,
                                                       gapLengthThreshold,
                                                       nThreshold,
                                                       atomic);
  rea->identifyFromLeftBreakpoint(top);
  return RearrangementPtr(rea);
}

GappedTopSegmentIteratorConstPtr MMapGenome::getGappedTopSegmentIterator(
  hal_index_t i, hal_size_t gapThreshold, bool atomic) const
{
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(i);  
  DefaultGappedTopSegmentIterator* gt = 
     new DefaultGappedTopSegmentIterator(top, gapThreshold, atomic);
  return GappedTopSegmentIteratorConstPtr(gt);
}

GappedBottomSegmentIteratorConstPtr MMapGenome::getGappedBottomSegmentIterator(
  hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
  bool atomic) const
{
  BottomSegmentIteratorConstPtr bot = getBottomSegmentIterator(i);  
  DefaultGappedBottomSegmentIterator* gb = 
     new DefaultGappedBottomSegmentIterator(bot, childIdx, gapThreshold, 
                                            atomic);
  return GappedBottomSegmentIteratorConstPtr(gb);
}

void MMapGenome::deleteSequenceCache()
{
  if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0)
  {
    map<hal_size_t, MMapSequence*>::iterator i;
    for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i)
    {
      delete i->second;
    }
    vector<MMapSequence*>::iterator z;
    for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z)
    {
      delete *z;
    }
  }
  else if (_sequenceNameCache.size() > 0)
  {
    map<string, MMapSequence*>::iterator i;
    for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i)
    {
      delete i->second;
    }
  }
  _sequencePosCache.clear();
  _zeroLenPosCache.clear();
  _sequenceNameCache.clear(); // I share my pointers with above. 
}

void MMapGenome::loadSequenceNameCache() const
{
  if (_sequenceNameCache.size() > 0)
  {
    return;
  }
  hal_size_t numSequences = _data->_numSequences;
  for (hal_size_t i = 0; i < numSequences; ++i)
  {
      MMapSequence* seq = 
          new MMapSequence(const_cast<MMapGenome*>(this),
                           i);

      _sequenceNameCache.insert(
          pair<string, MMapSequence*>(seq->getName(), seq));
  }
}


void MMapGenome::loadSequencePosCache() const
{
  if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0)
  {
    return;
  }
  hal_size_t totalReadLen = 0;
  loadSequenceNameCache();
  map<std::string, MMapSequence*>::iterator i;
  for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i)
  {
      if (i->second->getSequenceLength() > 0)
      {
          _sequencePosCache.insert(pair<hal_size_t, MMapSequence*>(
                                       i->second->getStartPosition() +
                                       i->second->getSequenceLength(), i->second));
          totalReadLen += i->second->getSequenceLength();
      }
      else
      {
          _zeroLenPosCache.push_back(i->second);
      }
  }
}

