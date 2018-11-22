rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/maf

libHalMaf_srcs = impl/halMafBed.cpp impl/halMafBlock.cpp impl/halMafExport.cpp \
    impl/halMafScanDimensions.cpp impl/halMafScanner.cpp impl/halMafScanReference.cpp \
    impl/halMafWriteGenomes.cpp
libHalMaf_objs = ${libHalMaf_srcs:%.cpp=${modObjDir}/%.o}
hal2maf_srcs = impl/hal2maf.cpp
hal2maf_objs = ${hal2maf_srcs:%.cpp=${modObjDir}/%.o}
maf2hal_srcs = impl/maf2hal.cpp
maf2hal_objs = ${maf2hal_srcs:%.cpp=${modObjDir}/%.o}
halMafTests_srcs = tests/halMafTests.cpp tests/halMafBlockTest.cpp tests/halMafExportTest.cpp
halMafTests_objs = ${halMafTests_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalMaf_srcs} ${hal2maf_srcs} ${maf2hal_srcs} ${halMafTests_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/hal2maf ${binDir}/maf2hal ${binDir}/halMafTests
otherLibs = ${libHalMaf} ${libHalLiftover} ${halApiTestSupportLibs}
inclSpec += -I${rootDir}/liftover/inc -I${halApiTestIncl}

all: libs progs
libs: ${libHalMaf}
progs: ${progs}

clean : 
	rm -f ${libHalMaf} ${objs} ${progs} ${depends} output

test: halMafTests hal2MafSmallMMapTest hal2MafSmallHdf5Test naiveLiftUpTests

halMafTests:
	${binDir}/halMafTests

naiveLiftUpTests:
	python2 -m pytest impl/naiveLiftUp.py

hal2MafSmallMMapTest:
	@mkdir -p output
	../bin/halRandGen --preset small --seed 0 --testRand --format mmap output/small.mmap.hal
	../bin/hal2maf output/small.mmap.hal output/small.mmap.maf
	diff tests/expected/small.mmap.maf output/small.mmap.maf

hal2MafSmallHdf5Test:
	@mkdir -p output
	../bin/halRandGen --preset small --seed 0 --testRand --format hdf5 output/small.hdf5.hal
	../bin/hal2maf output/small.hdf5.hal output/small.hdf5.maf
	diff tests/expected/small.hdf5.maf output/small.hdf5.maf

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

