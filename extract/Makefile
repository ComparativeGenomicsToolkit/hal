rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/extract

halExtract_srcs = impl/halExtract.cpp
halExtract_objs = ${halExtract_srcs:%.cpp=${modObjDir}/%.o}
halAlignedExtract_srcs = impl/halAlignedExtract.cpp
halAlignedExtract_objs = ${halAlignedExtract_srcs:%.cpp=${modObjDir}/%.o}
halMaskExtract_srcs = impl/halMaskExtractMain.cpp impl/halMaskExtractor.cpp
halMaskExtract_objs = ${halMaskExtract_srcs:%.cpp=${modObjDir}/%.o}
hal4dExtract_srcs = impl/hal4dExtractMain.cpp impl/hal4dExtract.cpp
hal4dExtract_objs = ${hal4dExtract_srcs:%.cpp=${modObjDir}/%.o}
halSingleCopyRegionsExtract_srcs = impl/halSingleCopyRegionsExtract.cpp
halSingleCopyRegionsExtract_objs = ${halSingleCopyRegionsExtract_srcs:%.cpp=${modObjDir}/%.o}
hal4dExtractTest_srcs = tests/hal4dExtractTest.cpp
hal4dExtractTest_objs = ${hal4dExtractTest_srcs:%.cpp=${modObjDir}/%.o} ${modObjDir}/impl/hal4dExtract.o
srcs = ${halExtract_srcs} ${halAlignedExtract_srcs} ${halMaskExtract_srcs} \
    ${hal4dExtract_srcs} ${halSingleCopyRegionsExtract_srcs} ${hal4dExtractTest_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halStats ${binDir}/halCoverage
inclSpec += -I${rootDir}/liftover/inc -I${halApiTestIncl}
otherLibs += ${halApiTestSupportLibs} ${libHalLiftover}
progs = ${binDir}/halExtract ${binDir}/halAlignedExtract ${binDir}/halMaskExtract \
    ${binDir}/hal4dExtract ${binDir}/halSingleCopyRegionsExtract ${binDir}/hal4dExtractTest

all: progs
libs:
progs: ${progs}

clean : 
	rm -f ${objs} ${progs} ${depends}
test: hal4dExtractTest
hal4dExtractTest:
	${binDir}/hal4dExtractTest 

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

