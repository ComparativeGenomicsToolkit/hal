/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <string>
#include <iostream>
#include "halAlignmentTest.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignment.h"
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

TempCreateAlignment::TempCreateAlignment(AlignmentPtr alignment) :
  _alignment(alignment)
{
  _path = getTempFile();
  _alignment->createNew(_path);
}

TempCreateAlignment::~TempCreateAlignment()
{
  _alignment->close();
  removeTempFile(const_cast<char*>(_path.c_str()));
}

TempReadAlignment::TempReadAlignment(AlignmentPtr alignment, 
                                     const string& path)
  : _path(path)
{
  alignment->open(_path, true);
  _alignment = alignment;
}

TempReadAlignment::~TempReadAlignment()
{
  _alignment->close();
  removeTempFile(const_cast<char*>(_path.c_str()));
}

void AlignmentTest::check()
{
  vector<AlignmentPtr> createInstances = getTestAlignmentInstances();
  vector<AlignmentPtr> readInstances = getTestAlignmentInstances();

  for (size_t i = 0; i < createInstances.size(); ++i)
  {
    TempCreateAlignment creater(createInstances[i]);
    createCallBack(creater._alignment);
    
    TempReadAlignment checker(readInstances[i], creater._path);
    checkCallBack(checker._alignment);
  }
}

