rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/chain

libHalChain_srcs = impl/halBlockViz.cpp impl/halChain.cpp
libHalChain_objs = ${libHalChain_srcs:%.cpp=${modObjDir}/%.o}
hal2chain_srcs = impl/hal2chain.cpp
hal2chain_objs = ${hal2chain_srcs:%.cpp=${modObjDir}/%.o}
halChainTests_srcs = tests/halChainTests.cpp tests/halChainGetBlocksTest.cpp
halChainTests_objs = ${halChainTests_srcs:%.cpp=${modObjDir}/%.o}
blockVizBed_srcs = tests/blockVizBed.cpp
blockVizBed_objs = ${blockVizBed_srcs:%.cpp=${modObjDir}/%.o}
blockVizMaf_srcs = tests/blockVizMaf.cpp
blockVizMaf_objs = ${blockVizMaf_srcs:%.cpp=${modObjDir}/%.o}
blockVizTest_srcs = tests/blockVizTest.cpp
blockVizTest_objs = ${blockVizTest_srcs:%.cpp=${modObjDir}/%.o}
blockVizTime_srcs = tests/blockVizTime.cpp
blockVizTime_objs = ${blockVizTime_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalChain_srcs} ${hal2chain_srcs} ${halChainTests_srcs} ${blockVizBed_srcs} \
    ${blockVizMaf_srcs} ${blockVizTest_srcs} ${blockVizTime_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
inclSpec += -I${rootDir}/liftover/inc -I${rootDir}/lod/inc -I${rootDir}/maf/inc -I${halApiTestIncl}
otherLibs += ${halApiTestSupportLibs} ${libHalChain} ${libHalLiftover} ${libHalLod} ${libHalMaf}
progs = ${binDir}/hal2chain ${binDir}/halChainTests ${binDir}/blockVizBed \
    ${binDir}/blockVizMaf ${binDir}/blockVizTest ${binDir}/blockVizTime \
    ${binDir}/halChainTests

testTmpDir = tests.tmp

all: libs progs
libs: ${libHalChain}
progs: ${progs}

clean: 
	rm -f ${libHalChain} ${objs} ${progs} ${depends} ${testTmpDir}

testHdf5Hal = ${testTmpDir}/small.haf5.hal
testMmapHal = ${testTmpDir}/small.mmap.hal

test: halChainTests blockVizHdf5Tests blockVizMmapTests

halChainTests:
	${binDir}/halChainTests

blockVizHdf5Tests: ${testHdf5Hal}
	${binDir}/blockVizTest ${testHdf5Hal} Genome_4 Genome_0 Genome_0_seq 0 3000


blockVizMmapTests: ${testMmapHal}
	${binDir}/blockVizTest ${testMmapHal} Genome_4 Genome_0 Genome_0_seq 0 3000

${testHdf5Hal}:
	@mkdir -p $(dir $@)
	${binDir}/halRandGen --preset small --seed 0 --format hdf5 $@

${testMmapHal}:
	@mkdir -p $(dir $@)
	${binDir}/halRandGen --preset small --seed 0 --format mmap $@

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

