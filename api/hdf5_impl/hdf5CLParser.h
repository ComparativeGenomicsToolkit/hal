/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5CLPARSER_H
#define _HDF5CLPARSER_H

#include <H5Cpp.h>
#include "halCLParser.h"
#include "halCLParserInstance.h"

namespace hal {

/** 
 * HDF5 extension of hal::CLParser
 */
class HDF5CLParser : public CLParser
{
public:
   ~HDF5CLParser();

   void applyToDCProps(H5::DSetCreatPropList& dcprops) const;
   void applyToAProps(H5::FileAccPropList& aprops) const;
   bool getInMemory() const;

   static const hsize_t DefaultChunkSize;
   static const hsize_t DefaultDeflate;
   static const hsize_t DefaultCacheMDCElems;
   static const hsize_t DefaultCacheRDCElems;
   static const hsize_t DefaultCacheRDCBytes;
   static const double DefaultCacheW0;
   static const bool DefaultInMemory;

protected:
   // Nobody creates this class except through the interface. 
   friend CLParserPtr hdf5CLParserInstance(bool createOptions);
   friend class HDF5Alignment;

   HDF5CLParser(bool createOptions);
};

}
#endif

