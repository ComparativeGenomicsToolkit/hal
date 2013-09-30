rootPath = ../
include ../include.mk

phastCflags = -I../../phast/include -I../../clapack/INCLUDE
phastLinkflags = ../../phast/lib/libphast.a ../../phast/lib/liblapack.a ../../phast/lib/libblaswr.a ../../clapack/F2CLIBS/libf2c.a
libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halReplaceGenome ${binPath}/halAppendSubtree ${binPath}/findDuplications ${binPath}/findRegionsExclusivelyInGroup

clean : 
	rm -f ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halReplaceGenome

${binPath}/halRemoveGenome: halRemoveGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halRemoveGenome halRemoveGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halAddToBranch: halAddToBranch.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAddToBranch halAddToBranch.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halReplaceGenome: halReplaceGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halReplaceGenome halReplaceGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halAppendSubtree: halAppendSubtree.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAppendSubtree halAppendSubtree.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

# TMP
${binPath}/findDuplications: findDuplications.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/findDuplications findDuplications.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/findRegionsExclusivelyInGroup: findRegionsExclusivelyInGroup.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/findRegionsExclusivelyInGroup findRegionsExclusivelyInGroup.cpp ${libPath}/halLib.a ${libPath}/halLiftover.a ${basicLibs}

${binPath}/ancestorsML: ancestorsML.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} ancestorsML.cpp ${basicLibs} ${libPath}/halLib.a ${phastLinkflags}

${binPath}/adjacenciesParsimony: adjacenciesParsimony.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/adjacenciesParsimony adjacenciesParsimony.cpp ${libPath}/halLib.a ${basicLibs}

markAncestors.o: markAncestors.cpp
	${cpp} ${cppflags} -I ${libPath} -c -o markAncestors.o markAncestors.cpp
