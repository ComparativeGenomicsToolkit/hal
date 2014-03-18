#ifndef __ANCESTORSMLBED_H_
#define __ANCESTORSMLBED_H_
#include <fstream>
#include <iostream>
#include "hal.h"
extern "C" {
#include "tree_model.h"
}
#include "ancestorsML.h"
#include "ancestorsMLBed.h"

using namespace hal;
using namespace std;

void AncestorsMLBed::visitLine()
{
  hal_index_t startPos, endPos;
  string sequenceName = _bedLine._chrName;
  startPos = _bedLine._start;
  endPos = _bedLine._end;
  Sequence *sequence = _genome->getSequence(sequenceName);
  if (sequence == NULL) {
    throw hal_exception("Sequence name not found!");
  }
  startPos += sequence->getStartPosition();
  endPos += sequence->getStartPosition();

  reEstimate(_mod, _alignment, _genome, startPos, endPos, _nameToId, _threshold, _writeHal, _printWrites);
}

#endif
