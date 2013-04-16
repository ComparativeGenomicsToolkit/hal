/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODMANAGER_H
#define _HALLODMANAGER_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"

namespace hal {

/** This is a container that keeps track of LOD alignments as generated
 * by halLodExtract.py
 */
class LodManager
{
public:
   
   LodManager();
   virtual ~LodManager();
   
   /** Load series of alignments specified in the lodPath file.  Options
    * from the given CLParser are applied if specified */
   void loadLODFile(const std::string& lodPath,
                    CLParserConstPtr options = CLParserConstPtr());

   /** Just use the given HAL file for everything.  Same as if we gave a
    * lodFile containing only "0 halPath"*/
   void loadSingeHALFile(const std::string& halPath,
                         CLParserConstPtr options = CLParserConstPtr());

   AlignmentConstPtr getAlignment(hal_size_t queryLength) const;
   
protected:

   void checkMap(const std::string& lodPath) const;
   typedef std::map<hal_size_t, AlignmentConstPtr> AlignmentMap;

   AlignmentMap _map;

};

HAL_FORWARD_DEC_CLASS(LodManager)


}

#endif
