rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/liftover

libHalLiftover_srcs = impl/halBedLine.cpp impl/halBedScanner.cpp impl/halBlockLiftover.cpp \
    impl/halBlockMapper.cpp impl/halColumnLiftover.cpp impl/halLiftover.cpp \
    impl/halWiggleLiftover.cpp impl/halWiggleLoader.cpp impl/halWiggleScanner.cpp
libHalLiftover_objs = ${libHalLiftover_srcs:%.cpp=${modObjDir}/%.o}
halLiftover_srcs = impl/halLiftoverMain.cpp
halLiftover_objs = ${halLiftover_srcs:%.cpp=${modObjDir}/%.o}
halWiggleLiftover_srcs = impl/halWiggleLiftoverMain.cpp
halWiggleLiftover_objs = ${halWiggleLiftover_srcs:%.cpp=${modObjDir}/%.o}
halLiftoverTests_srcs = tests/halLiftoverTests.cpp
halLiftoverTests_objs = ${halLiftoverTests_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalLiftover_srcs} ${halLiftover_srcs} ${halWiggleLiftover_srcs} ${halLiftover_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halLiftover ${binDir}/halWiggleLiftover ${binDir}/halLiftoverTests
otherLibs += ${libHalLiftover} ${halApiTestSupportLibs}

# tests use api/tests/halAlignmentTest
inclSpec += -I${rootDir}/api/tests
otherLibs += ${halApiTestSupportLib}

all : libs progs
libs: ${libHalLiftover}
progs: ${progs}

clean: 
	rm -rf ${libHalLiftover} ${objs} ${progs} ${depends} output

test: unitTests halLiftoverBed12Test halLiftoverPsl12Test \
	 halLiftoverBed3Test halLiftoverPsl3Test

unitTests:
	${binDir}/halLiftoverTests 

halLiftoverBed12Test: output/small.hdf5.hal
	${binDir}/halLiftover output/small.hdf5.hal Genome_0 tests/input/test1.bed12 Genome_2 output/$@.bed
	diff -u tests/expected/$@.bed output/$@.bed

halLiftoverPsl12Test: output/small.hdf5.hal
	${binDir}/halLiftover --outPSL output/small.hdf5.hal Genome_0 tests/input/test1.bed12 Genome_2 output/$@.psl
	diff -u tests/expected/$@.psl output/$@.psl

halLiftoverBed3Test: output/small.hdf5.hal
	${binDir}/halLiftover output/small.hdf5.hal Genome_0 tests/input/test1.bed3 Genome_2 output/$@.bed
	diff -u tests/expected/$@.bed output/$@.bed

halLiftoverPsl3Test: output/small.hdf5.hal
	${binDir}/halLiftover --outPSL output/small.hdf5.hal Genome_0 tests/input/test1.bed3 Genome_2 output/$@.psl
	diff -u tests/expected/$@.psl output/$@.psl

output/small.hdf5.hal: ../bin/halRandGen
	@mkdir -p output
	../bin/halRandGen --preset small --seed 0 --testRand --format hdf5 output/small.hdf5.hal

../bin/halRandGen:
	cd ../randgen && ${MAKE}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

