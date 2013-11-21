rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/hal2chain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTestsAll = $(wildcard test/*.cpp) 
libTests=$(subst test/blockInterpolateTest.cpp,,${libTestsAll})
libTestsHeaders = $(wildcard test/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${libPath}/halChain.a ${binPath}/hal2chain test/blockVizTest test/blockVizBed test/blockVizTime test/blockVizMaf ${binPath}/halChainTests 

clean : 
	rm -f ${libPath}/halChain.a ${libPath}/*.h ${binPath}/hal2chain test/blockVizTest test/BlockVizBed test/blockVizTime test/blockVizMaf 

${libPath}/halChain.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halChain.a *.o
	ranlib halChain.a 
	rm *.o
	mv halChain.a ${libPath}/

${binPath}/hal2chain : impl/hal2chain.cpp ${libPath}/halMaf.a ${libPath}/halChain.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/hal2chain impl/hal2chain.cpp ${libPath}/halMaf.a ${libPath}/halChain.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

test/blockVizTest: test/blockVizTest.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockVizTest test/blockVizTest.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

test/blockVizMaf: test/blockVizMaf.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockVizMaf test/blockVizMaf.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

test/blockVizBed: test/blockVizBed.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockVizBed test/blockVizBed.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

test/blockVizTime: test/blockVizTime.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockVizTime test/blockVizTime.c ${libPath}/halChain.a ${libPath}/halLod.a ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halChainTests : ${libTests} ${libTestsHeaders} ${libTestsCommon} ${libTestsHeadersCommon} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -I ../api/tests -o ${binPath}/halChainTests ${libHalTests} ${libTests} ${libPath}/halLib.a ${libPath}/halMaf.a ${libPath}/halChain.a ${libPath}/halLiftover.a ${basicLibs}

