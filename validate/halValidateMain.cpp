/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halStats.h"

using namespace std;
using namespace hal;

static string checkOptions(int argc, char** argv)
{
  if (argc != 2)
  {
    cerr << "Usage: halValidate <path of hal alignment file>" << endl;
    exit(1);
  }
  return argv[1];
}

int main(int argc, char** argv)
{
  string path = checkOptions(argc, argv);
  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(path);
    validateAlignment(alignment);
  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  cout << "\nFile valid" << endl;
  
  return 0;
}
