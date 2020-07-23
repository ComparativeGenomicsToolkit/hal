/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODMANAGER_H
#define _HALLODMANAGER_H

#include "hal.h"
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace hal {
    class CLParser;

    /** This is a container that keeps track of LOD alignments as generated
     * by halLodExtract.py
     */
    class LodManager {
      public:
        LodManager();
        virtual ~LodManager();

        /** Load series of alignments specified in the lodPath file.  Options
         * from the given CLParser are applied if specified.
         *
         * If the paths of the HAL files are relative (do not begin with /) then
         * they will be concatenated to the directory of lodPath.  If they
         * are absolute (beginning with /) then they will be opened directly.
         * Paths that contain ":/" are assumed to be
         * web addressed of some sort and considered absolute. */
        void loadLODFile(const std::string &lodPath, const CLParser *options = NULL);

        /** Just use the given HAL file for everything.  Same as if we gave a
         * lodFile containing only "0 halPath"*/
        void loadSingeHALFile(const std::string &halPath, const CLParser *options = NULL);

        AlignmentConstPtr getAlignment(hal_size_t queryLength, bool needDNA);

        /** Check if query length corresponds to LOD 0 (ie original HAL) */
        bool isLod0(hal_size_t queryLenth) const;

        /** Any query greater than this is disabled */
        hal_size_t getMaxQueryLength() const;

        /** Maximum age of a URL in seconds such that we dont try to
         * preload headers for all the HAL files */
        static const unsigned long MaxAgeSec;

        /** Token that specifies upper limit for LODs, that sits in path field */
        static const std::string MaxLodToken;

      private:
        std::string resolvePath(const std::string &lodPath, const std::string &halPath);
        void checkMap(const std::string &lodPath);
        void checkAlignment(hal_size_t minQuery, const std::string &path, AlignmentConstPtr alignment);
        void preloadAlignments();

        typedef std::pair<std::string, AlignmentConstPtr> PathAlign;
        typedef std::map<hal_size_t, PathAlign> AlignmentMap;

        const CLParser *_options;
        AlignmentMap _map;
        hal_size_t _maxLodLowerBound;
    };

    inline hal_size_t LodManager::getMaxQueryLength() const {
        return _maxLodLowerBound - 1;
    }

    HAL_FORWARD_DEC_CLASS(LodManager)
}

#endif
// Local Variables:
// mode: c++
// End:
