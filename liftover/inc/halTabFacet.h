/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// based on:
//http://stackoverflow.com/questions/9085471/is-there-a-flag-to-make-istream-treat-only-tabs-as-delimiters

#include <locale>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>    

// This is my facet:
// It is designed to treat only <tab> as whitespace
class TabSepFacet: public std::ctype<char>
{
public:
   typedef std::ctype<char>   base;
   typedef base::char_type    char_type;

   TabSepFacet(std::locale const& l) : base(table)
   {
     // Get the ctype facet of the current locale
     std::ctype<char> const& defaultCType = 
        std::use_facet<std::ctype<char> >(l);

     // Copy the default flags for each character from the current facet
     static char data[256];
     for(int loop = 0; loop < 256; ++loop) 
     {
       data[loop] = loop;
     }
     defaultCType.is(data, data+256, table);
     
     // Remove the space character
     if (table[' '] & base::space)
     {   
       table[' '] ^= base::space;
     }
    
     // Add the tab character
     table['\t'] |= base::space;
   }
private:
   base::mask table[256];
};
