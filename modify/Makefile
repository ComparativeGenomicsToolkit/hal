rootPath = ../
include ../include.mk

libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halAddToNode ${binPath}/halReplaceGenome ${binPath}/halAppendSubtree

clean : 
	rm ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halAddToNode ${binPath}/halReplaceGenome

${binPath}/halRemoveGenome: halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halRemoveGenome halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halAddToBranch: halAddToBranch.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAddToBranch halAddToBranch.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halAddToNode: halAddToNode.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAddToNode halAddToNode.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halReplaceGenome: halReplaceGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halReplaceGenome halReplaceGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halAppendSubtree: halAppendSubtree.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAppendSubtree halAppendSubtree.cpp ${libPath}/halLib.a ${basicLibs}