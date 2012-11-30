/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <locale>

using namespace std;

static size_t countBlock(vector<string> block, size_t rows);
static char reverseComplement(char c);
static void reverseComplement(std::string& s);

// quick hacky tool to quickly count mutations in maf.  for benchmark baselines
// g++ -O3 mafMutations.cpp -o mafMutations
int main(int argc, char** argv)
{
  if (argc != 2)
  {
    cerr << "usage : mafMutations <mafFile>" << endl;
    return 1;
  }
  
  string mafPath = argv[1];
  ifstream mafFile(mafPath.c_str());
  
  if (!mafFile)
  {
    cerr << "problem reading " << mafPath << endl;
    return 1;
  }

  // copied from halMafScanner.cpp
  string buffer;
  string sequenceName;
  size_t startPosition, length, srcLength;
  char strand;
  vector<string> block;
  size_t count = 0;
  size_t rows = 0;
  while (!mafFile.eof() && mafFile.good())
  {
    buffer.clear();
    mafFile >> buffer;
    if (buffer == "a")
    {
      count += countBlock(block, rows);
      rows = 0;
    }
    else if (buffer == "s")
    {
      if (rows + 1 > block.size())
      {
        block.resize(rows + 1);
      }
      mafFile >> sequenceName >> startPosition >> length 
              >> strand >> srcLength >> block[rows];
      if (strand == '-')
      {
        reverseComplement(block[rows]);
      }
      if (!mafFile.good())
      {
        cerr << "problem parsing " << mafPath << endl;
        return 1;
      }
      ++rows;
    }
    else
    {
      while (!mafFile.eof() && !mafFile.bad() && mafFile.peek() != '\n')
      {
        mafFile.get();
      }
    }
  }
  count += countBlock(block, rows);
  cout << "count = " << count << endl;
}

size_t countBlock(vector<string> block, size_t rows)
{
  size_t count = 0;
  for (size_t i = 1; i < rows; ++i)
  {
    if (block[i].length() != block[0].length())
    {
      cerr << "block[" << i << "] = " << block[i] 
           << "\nblock[0] = " << block[0] << endl;
      cerr << "block parse error " << endl;
      exit(1);
    }
    for (size_t j = 0; j < block[i].length(); ++j)
    {
      char a = block[0][j];
      char b = block[i][j];
      if (a != '-' && b != '-')
      {
        a = toupper(a);
        b = toupper(b);
        if (a != b)
        {
          ++count;
        }
      }
    }
  }
  return count;
}

void reverseComplement(std::string& s)
{
  if (!s.empty())
  {
    size_t j = s.length() - 1;
    size_t i = 0;
    char buf;
    do
    {
      while (j > 0 && s[j] == '-')
      {
        --j;
      }
      while (i < s.length() - 1 && s[i] == '-')
      {
        ++i;
      }
      
      if (i >= j || s[i] == '-' || s[j] == '-')
      {
        if (i == j && s[i] != '-')
        {
          s[i] = reverseComplement(s[i]);
        }
        break;
      }

      buf = reverseComplement(s[i]);
      s[i] = reverseComplement(s[j]);
      s[j] = buf;

      ++i;
      --j;
    } while (true);
  }
}

char reverseComplement(char c)
{
  switch (c)
  {
  case 'A' : return 'T'; 
  case 'a' : return 't'; 
  case 'C' : return 'G'; 
  case 'c' : return 'g';
  case 'G' : return 'C';
  case 'g' : return 'c';
  case 'T' : return 'A';
  case 't' : return 'a';
  default : break;
  }
  return c;
}
