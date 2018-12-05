/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>

#include "halCLParser.h"
#include "hal.h"
#include "halRandomData.h"
#include "halRandNumberGen.h"

using namespace std;
using namespace hal;

struct RandOptions {
   double _meanDegree;
   double _maxBranchLength;
   hal_size_t _minGenomes;
   hal_size_t _maxGenomes;
   hal_size_t _minSegmentLength;
   hal_size_t _maxSegmentLength;
   hal_size_t _minSegments;
   hal_size_t _maxSegments;
   int _seed;
   bool _testRand;
   string _halFile;
};

static const RandOptions defaultSm =  {0.75,  0.1,  2,   5,  250,   1000,      5,     10, -1, false, ""};
static const RandOptions defaultMed = {1.25,  0.7,  8,  20,  500,   2000,    100,    500, -1, false, ""};
static const RandOptions defaultBig = {2.00,  0.7, 20,  50,  1000,  8000,    400,   5000, -1, false, ""};
static const RandOptions defaultLrg = {2.00,  1.0, 50, 100,  5000, 10000,  10000,  50000, -1, false, ""};

static void initParser(CLParser& optionsParser) {
    optionsParser.setDescription("Generate a random HAL alignment file");
    optionsParser.addOption("preset", "one of small, medium, big, large [medium]", "medium");
    optionsParser.addOption("meanDegree", "[" + std::to_string(defaultMed._meanDegree) + "]", defaultMed._meanDegree);
    optionsParser.addOption("maxBranchLength", "[" + std::to_string(defaultMed._maxBranchLength) + "]", defaultMed._maxBranchLength);
    optionsParser.addOption("minGenomes", "[" + std::to_string(defaultMed._minGenomes) + "]", defaultMed._minGenomes);
    optionsParser.addOption("maxGenomes", "[" + std::to_string(defaultMed._maxGenomes) + "]", defaultMed._maxGenomes);
    optionsParser.addOption("minSegmentLength", "[" + std::to_string(defaultMed._minSegmentLength) + "]", defaultMed._minSegmentLength);
    optionsParser.addOption("maxSegmentLength", "[" + std::to_string(defaultMed._maxSegmentLength) + "]", defaultMed._maxSegmentLength);
    optionsParser.addOption("maxSegments", "[" + std::to_string(defaultMed._maxSegments) + "]", defaultMed._maxSegments);
    optionsParser.addOption("minSegments", "[" + std::to_string(defaultMed._minSegments) + "]", defaultMed._minSegments);
    optionsParser.addOption("seed", "random number seed ", -1);
    optionsParser.addOptionFlag("testRand", "use portable random number generator", false);
    optionsParser.addArgument("halFile", "path to toutput HAL alignment file");
}

template <typename T>
void updateOption(const CLParser* parser,
                  const string& name,
                  T& val) {
    if (parser->specifiedOption(name)) {
        val = parser->getOption<T>(name);
    }
}

static RandOptions getPresetDefault(const CLParser* optionsParser) {
    const string preset = optionsParser->getOption<string>("preset");
    if (preset == "small") {
        return defaultSm;
    } else if (preset == "medium") {
        return defaultMed;
    } else if (preset == "big") {
        return defaultBig;
    } else if (preset == "large") {
        return defaultLrg;
    } else {
        throw hal_exception(" invalid --preset value: " + preset);
    }
}

static RandOptions parseProgOptions(const CLParser* optionsParser)
{
    RandOptions options = getPresetDefault(optionsParser);
    updateOption(optionsParser, "meanDegree", options._meanDegree);
    updateOption(optionsParser, "maxBranchLength", options._maxBranchLength);
    updateOption(optionsParser, "minGenomes", options._minGenomes);
    updateOption(optionsParser, "maxGenomes", options._maxGenomes);
    updateOption(optionsParser, "minSegmentLength", options._minSegmentLength);
    updateOption(optionsParser, "maxSegmentLength", options._maxSegmentLength);
    updateOption(optionsParser, "minSegments", options._minSegments);
    updateOption(optionsParser, "maxSegments", options._maxSegments);
    updateOption(optionsParser, "seed", options._seed);
    if (optionsParser->getFlag("testRand")) {
        options._testRand = true;
    }
    options._halFile = optionsParser->getArgument<string>("halFile");
    return options;
}
 
int main(int argc, char** argv)
{
    CLParser optionsParser(CREATE_ACCESS);
    initParser(optionsParser);
    RandOptions options;
    try {
        optionsParser.parseOptions(argc, argv);
        options = parseProgOptions(&optionsParser);
    } catch (hal_exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    RandNumberGen rng(false, options._seed);
        
  try
  {
      AlignmentPtr alignment(openHalAlignment(options._halFile, &optionsParser, hal::CREATE_ACCESS));
    // call the crappy unit-test simulator 
      createRandomAlignment(rng, alignment.get(), options._meanDegree,
                            options._maxBranchLength, options._minGenomes,
                            options._maxGenomes, options._minSegmentLength,
                            options._maxSegmentLength, options._minSegments,
                            options._maxSegments);
    
    alignment->close();
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
  
  return 0;
}
