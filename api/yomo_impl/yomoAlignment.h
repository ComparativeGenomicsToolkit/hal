/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOALIGNMENT_H
#define _YOMOALIGNMENT_H

#include "halAlignmentInstance.h"
#include "yomoAlignment.h"
#include "yomoGenome.h"
#include "yomoMetaData.h"
#include <H5Cpp.h>
#include <map>

typedef struct _stTree stTree;

namespace hal {

    class YomoGenome;
    /**
     * YOMO implementation of hal::Alignment
     */
    class YomoAlignment : public Alignment {
      public:
        /* check if first bit of file has YOMO header */
        static bool isYomoFile(const std::string &initialBytes);

        YomoAlignment(const std::string &alignmentPath, unsigned mode, const H5::FileCreatPropList &fileCreateProps,
                      const H5::FileAccPropList &fileAccessProps, const H5::DSetCreatPropList &datasetCreateProps,
                      bool inMemory = false);
        YomoAlignment(const std::string &alignmentPath, unsigned mode, const CLParser *parser);
        ~YomoAlignment();

        static void defineOptions(CLParser *parser, unsigned mode);

        void close();

        const std::string &getStorageFormat() const {
            return STORAGE_FORMAT_YOMO;
        }

        Genome *addLeafGenome(const std::string &name, const std::string &parentName, double branchLength);

        Genome *addRootGenome(const std::string &name, double branchLength);

        void removeGenome(const std::string &name);

        Genome *insertGenome(const std::string &name, const std::string &parentName, const std::string &childName,
                             double upperBranchLength);

        const Genome *openGenome(const std::string &name) const;

        Genome *openGenome(const std::string &name);

        void closeGenome(const Genome *genome) const;

        std::string getRootName() const;

        std::string getParentName(const std::string &name) const;

        void updateBranchLength(const std::string &parentName, const std::string &childName, double length);

        double getBranchLength(const std::string &parentName, const std::string &childName) const;

        std::vector<std::string> getChildNames(const std::string &name) const;

        std::vector<std::string> getLeafNamesBelow(const std::string &name) const;

        hal_size_t getNumGenomes() const;

        MetaData *getMetaData();

        const MetaData *getMetaData() const;

        std::string getNewickTree() const;

        std::string getVersion() const;

        bool isReadOnly() const;

        void replaceNewickTree(const std::string &newNewickString);

      private:
        // FIXME: should these be private?
        void loadTree();
        void writeTree();
        void writeVersion();
        void addGenomeToTree(const std::string &name, const std::pair<std::string, double> &parentName,
                             const std::vector<std::pair<std::string, double>> &childNames);

      private:
        YomoAlignment() {
        }
        YomoAlignment(const YomoAlignment &) {
        }
        void initializeFromOptions(const CLParser *parser);
        void create();
        void open();
        void setInMemory();

      public:
        static const hsize_t DefaultChunkSize;
        static const hsize_t DefaultCompression;
        static const hsize_t DefaultCacheMDCElems;
        static const hsize_t DefaultCacheRDCElems;
        static const hsize_t DefaultCacheRDCBytes;
        static const double DefaultCacheW0;
        static const bool DefaultInMemory;

        static const H5std_string MetaGroupName;
        static const H5std_string TreeGroupName;
        static const H5std_string GenomesGroupName;
        static const H5std_string VersionGroupName;
        const std::string _alignmentPath;

      private:
        unsigned _mode;

        H5::H5File *_file;
        int _flags;
        bool _inMemory;
        H5::FileCreatPropList _cprops;
        H5::FileAccPropList _aprops;
        H5::DSetCreatPropList _dcprops;
        YOMOMetaData *_metaData;
        stTree *_tree;
        mutable std::map<std::string, stTree *> _nodeMap;
        bool _dirty;
        mutable std::map<std::string, YomoGenome *> _openGenomes;
    };
}
#endif

// Local Variables:
// mode: c++
// End:
