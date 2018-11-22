#include "mmapBottomSegment.h"
#include "mmapTopSegment.h"
#include "halDnaIterator.h"

using namespace hal;

void MMapBottomSegment::setCoordinates(hal_index_t startPos, hal_size_t length)
{
  if (_genome && (startPos >= (hal_index_t)_genome->getSequenceLength() || 
                  startPos + length > _genome->getSequenceLength()))
  {
    throw hal_exception("Trying to set top segment coordinate out of range");
  }

  _data->setStartPosition(startPos);
  getNextData()->setStartPosition(startPos + length);
}

hal_offset_t MMapBottomSegment::getTopParseOffset() const
{
  assert(_index >= 0);
  hal_offset_t offset = 0;
  hal_index_t topIndex = getTopParseIndex();
  if (topIndex != NULL_INDEX)
  {
      MMapTopSegment ts(getMMapGenome(), topIndex);
    assert(ts.getStartPosition() <= getStartPosition());
    assert((hal_index_t)(ts.getStartPosition() + ts.getLength()) 
           >= getStartPosition());
    offset = getStartPosition() - ts.getStartPosition();
  }
  return offset;
}

void MMapBottomSegment::print(std::ostream &os) const {
    os << "MMapBottomSegment" << getStartPosition() << " " << getEndPosition() << std::endl;
}
