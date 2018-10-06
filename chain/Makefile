rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/chain

libHalChain_srcs = impl/halBlockViz.cpp impl/halChain.cpp
libHalChain_objs = ${libHalChain_srcs:%.cpp=${modObjDir}/%.o}
hal2chain_srcs = impl/hal2chain.cpp
hal2chain_objs = ${hal2chain_srcs:%.cpp=${modObjDir}/%.o}
halChainTests_srcs = test/halChainTests.cpp test/halChainGetBlocksTest.cpp
halChainTests_objs = ${halChainTests_srcs:%.cpp=${modObjDir}/%.o}
blockVizBed_srcs = test/blockVizBed.cpp
blockVizBed_objs = ${blockVizBed_srcs:%.cpp=${modObjDir}/%.o}
blockVizMaf_srcs = test/blockVizMaf.cpp
blockVizMaf_objs = ${blockVizMaf_srcs:%.cpp=${modObjDir}/%.o}
blockVizTest_srcs = test/blockVizTest.cpp
blockVizTest_objs = ${blockVizTest_srcs:%.cpp=${modObjDir}/%.o}
blockVizTime_srcs = test/blockVizTime.cpp
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

all: libs progs
libs: ${libHalChain}
progs: ${progs}

clean: 
	rm -f ${libHalChain} ${objs} ${progs} ${depends}
test: halChainTests

halChainTests:
	${binDir}/halChainTests 

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

