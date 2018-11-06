#ifndef _MMAPDNAITERATOR_H
#define _MMAPDNAITERATOR_H
#include "halDNAIterator.h"
#include "halCommon.h"
namespace hal {
class MMapDNAIterator : public DNAIterator
{
public:
    MMapDNAIterator(MMapGenome* genome, hal_index_t index) :
        _index(index),
        _genome(genome),
        _reversed(false) {
    }

    char getChar() const;
    void setChar(char c);
    void toLeft() const;
    void toRight() const;
    void jumpTo(hal_size_t index) const;
    void toReverse() const;
    bool getReversed() const;
    void setReversed(bool reversed) const;
    const Genome* getGenome() const;
    Genome* getGenome();
    const Sequence* getSequence() const;
    hal_index_t getArrayIndex() const;

    bool equals(DNAIteratorPtr& other) const;
    bool leftOf(DNAIteratorPtr& other) const;

    void readString(std::string& outString, hal_size_t length) const;

    void writeString(const std::string& inString, hal_size_t length);

    bool inRange() const;

private:
    mutable hal_index_t _index;
    mutable MMapGenome *_genome;
    mutable bool _reversed;
};

inline bool MMapDNAIterator::inRange() const
{
  return _index >= 0 && 
     _index < (hal_index_t)_genome->getSequenceLength();
}

inline char MMapDNAIterator::getChar() const
{
  assert(inRange() == true);
  char c = *_genome->getDNA(_index, 1);
  if (_reversed)
  {
    c = reverseComplement(c);
  }
  return c;
}

inline void MMapDNAIterator::setChar(char c)
{
  if (inRange() == false) 
  {
    throw hal_exception("Trying to set character out of range");
  }
  else if (isNucleotide(c) == false)
  {
    throw hal_exception(std::string("Trying to set invalid character: ") + c);
  }
  if (_reversed)
  {
    c = reverseComplement(c);
  }
  *_genome->getDNA(_index, 1) = c;
}

inline void MMapDNAIterator::toLeft() const
{
  _reversed ? ++_index : --_index;
}

inline void MMapDNAIterator::toRight() const
{
  _reversed ? --_index : ++_index;
}

inline void MMapDNAIterator::jumpTo(hal_size_t index) const
{
  _index = index;
}

inline void MMapDNAIterator::toReverse() const
{
  _reversed = !_reversed;
}

inline bool MMapDNAIterator::getReversed() const
{
  return _reversed;
}

inline void MMapDNAIterator::setReversed(bool reversed) const
{
  _reversed = reversed;
}

inline const Genome* MMapDNAIterator::getGenome() const
{
  return _genome;
}

inline Genome* MMapDNAIterator::getGenome()
{
  return _genome;
}

inline const Sequence* MMapDNAIterator::getSequence() const
{
  return _genome->getSequenceBySite(_index);
}

inline hal_index_t MMapDNAIterator::getArrayIndex() const
{
  return _index;
}

inline bool MMapDNAIterator::equals(DNAIteratorPtr& other) const
{
  const MMapDNAIterator* mmOther = reinterpret_cast<
     const MMapDNAIterator*>(other.get());
  assert(_genome == mmOther->_genome);
  return _index == mmOther->_index;
}

inline bool MMapDNAIterator::leftOf(DNAIteratorPtr& other) const
{
  const MMapDNAIterator* mmOther = reinterpret_cast<
     const MMapDNAIterator*>(other.get());
  assert(_genome == mmOther->_genome);
  return _index < mmOther->_index;
}

inline void MMapDNAIterator::readString(std::string& outString,
                                        hal_size_t length) const
{
  assert(length == 0 || inRange() == true);
  outString.resize(length);

  for (hal_size_t i = 0; i < length; ++i)
  {
    outString[i] = getChar();
    toRight();
  }
}

inline void MMapDNAIterator::writeString(const std::string& inString,
                                         hal_size_t length)
{
  assert(length == 0 || inRange() == true);
  for (hal_size_t i = 0; i < length; ++i)
  {
    setChar(inString[i]);
    toRight();
  }
}

}
#endif
