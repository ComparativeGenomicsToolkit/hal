
#ifndef _RANDNUMBERGEN_H
#define _RANDNUMBERGEN_H

#include <ctime>
#include <math.h>
#include <random>
#include <stdlib.h>

namespace hal {

    /**
     * Interface to random number generator.  Normally, this is just a
     * pass-through interface to C++ random facility.  However, this can run in
     * test mode, using an internal random number generator that will produce the
     * same results on all systems.  This is used for tests so that they can be
     * reproduced on any system.  Internal random number code originally written
     * by William Stafford Noble.
     */
    class RandNumberGen {
      private:
        /* Maximum value using in generating test random number. */
        static const int TEST_RAND_MAX = 4096;

        /* Are we running in test mode? */
        bool _testMode;

        /* running seed for test mode */
        int _seed;

        std::mt19937 rng;

        /*
         * Get a pseudo-random number reproducable on all systems
         */
        double getTestRand() {
            _seed = abs((_seed / 3) * _seed + 7718);
            return double(_seed % TEST_RAND_MAX) / double(TEST_RAND_MAX);
        }

      public:
        /**
         * Constructor, if seed >=0 is used to set the seed.
         */
        RandNumberGen(bool testMode = false, int seed = -1) : _testMode(testMode), _seed(0) {
            if (_testMode) {
                if (_seed >= 0) {
                    _seed = seed;
                } else {
                    _seed = time(NULL);
                }
            } else {
                if (_seed >= 0) {
                    rng.seed(seed);
                }
            }
        }

        /**
         * Get a random number from 0.0 and 1.0.
         */
        inline double getRand() {
            if (_testMode) {
                return getTestRand();
            } else {
                std::uniform_real_distribution<double> dist;
                return dist(rng);
            }
        }

        /**
         * Get a random positive double in the specified range.
         */
        inline int getRandDouble(double minVal, double maxVal) {
            if (maxVal < minVal) {
                // This is useful for generating `at least' lengths,
                // even if them max is below this.
                maxVal = minVal;
            }
            return (getRand() * (maxVal - minVal)) + minVal;
        }

        /**
         * Get a random positive integer in the specified range.
         */
        inline int getRandInt(int minVal, int maxVal) {
            if (maxVal < minVal) {
                // This is useful for generating `at least' lengths,
                // even if them max is below this.
                maxVal = minVal;
            }
            double rnum = getRand() * double(maxVal - minVal);
            if ((rnum - floor(rnum)) >= 0.5) {
                return minVal + int(ceil(rnum));
            } else {
                return minVal + int(floor(rnum));
            }
        }

        /**
         * Get a random zero or positive integer less or equal to the specified
         * max.
         */
        inline int getRandInt(int maxVal) {
            return getRandInt(0, maxVal);
        }
    };
}
#endif

/*
 * Local Variables:
 * mode: c++
 * End:
 */
