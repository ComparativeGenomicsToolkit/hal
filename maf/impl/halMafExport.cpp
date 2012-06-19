/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halMafExport.h"

using namespace std;
using namespace hal;

MafExport::MafExport()
{

}

MafExport::~MafExport()
{

}

void MafExport::convert(AlignmentConstPtr alignment,
                        ostream& mafStream)
{
  _alignment = alignment;
  _mafStream = &mafStream;
}
