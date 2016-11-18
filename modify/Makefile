rootPath = ../
include ../include.mk

phastCflags = -I../../phast/include -I../../clapack/INCLUDE -I../../phast/src/lib/pcre
phastLinkflags = ../../phast/lib/libphast.a ../../phast/lib/liblapack.a ../../phast/lib/libblaswr.a ../../clapack/F2CLIBS/libf2c.a
libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})
targets = ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halReplaceGenome ${binPath}/halAppendSubtree ${binPath}/findRegionsExclusivelyInGroup ${binPath}/halUpdateBranchLengths ${binPath}/halWriteNucleotides ${binPath}/halSetMetadata ${binPath}/halRenameGenomes
phastTargets = ${binPath}/ancestorsML
ifdef ENABLE_PHYLOP
all : ${targets} ${phastTargets}
else
all : ${targets}
endif

clean : 
	rm -f ${targets}

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
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -c ancestorsMLBed.cpp -o ancestorsMLBed.o ${basicLibs} ${libPath}/halLib.a ${phastLinkflags}
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -o ${binPath}/ancestorsML ancestorsML.cpp ${basicLibs} ancestorsMLBed.o ${libPath}/halLib.a ${libPath}/halLiftover.a ${phastLinkflags}

${binPath}/adjacenciesParsimony: adjacenciesParsimony.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/adjacenciesParsimony adjacenciesParsimony.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halUpdateBranchLengths: halUpdateBranchLengths.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halUpdateBranchLengths halUpdateBranchLengths.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halWriteNucleotides: halWriteNucleotides.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halWriteNucleotides halWriteNucleotides.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halSetMetadata: halSetMetadata.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halSetMetadata halSetMetadata.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halRenameGenomes: halRenameGenomes.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o $@ $^

markAncestors.o: markAncestors.cpp
	${cpp} ${cppflags} -I ${libPath} -c -o markAncestors.o markAncestors.cpp
