#include "mmapBottomSegment.h"
#include "mmapTopSegment.h"
#include "mmapDNAIterator.h"

using namespace hal;

void MMapBottomSegment::setCoordinates(hal_index_t startPos, hal_size_t length)
{
  if (_genome && (startPos >= (hal_index_t)_genome->getSequenceLength() || 
                  startPos + length > _genome->getSequenceLength()))
  {
    throw hal_exception("Trying to set top segment coordinate out of range");
  }

  _data->setStartPosition(startPos);
  (_data + 1)->setStartPosition(startPos + length);
}

hal_offset_t MMapBottomSegment::getTopParseOffset() const
{
  assert(_index >= 0);
  hal_offset_t offset = 0;
  hal_index_t topIndex = getTopParseIndex();
  if (topIndex != NULL_INDEX)
  {
    MMapTopSegment ts(_genome, topIndex);
    assert(ts.getStartPosition() <= getStartPosition());
    assert((hal_index_t)(ts.getStartPosition() + ts.getLength()) 
           >= getStartPosition());
    offset = getStartPosition() - ts.getStartPosition();
  }
  return offset;
}

void MMapBottomSegment::getString(std::string& outString) const
{
  MMapDNAIterator di(const_cast<MMapGenome*>(_genome), getStartPosition());
  di.readString(outString, getLength());
}

bool MMapBottomSegment::isMissingData(double nThreshold) const
{
  if (nThreshold >= 1.0)
  {
    return false;
  }
  MMapDNAIterator di(const_cast<MMapGenome*>(_genome), getStartPosition());
  size_t length = getLength();
  size_t maxNs = nThreshold * (double)length;
  size_t Ns = 0;
  char c;
  for (size_t i = 0; i < length; ++i, di.toRight())
  {
    c = di.getChar();
    if (c == 'N' || c == 'n')
    {
      ++Ns;
    }
    if (Ns > maxNs)
    {
      return true;
    }
    if ((length - i) < (maxNs - Ns))
    {
      break;
    }
  }
  return false;
}
