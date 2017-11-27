rootPath = ../
include ../include.mk

phastCflags = ${phyloPcppflags}
phastLinkflags = ${phyloPlibs}
libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})
targets = ${binPath}/halRemoveGenome ${binPath}/halAddToBranch ${binPath}/halReplaceGenome ${binPath}/halAppendSubtree ${binPath}/findRegionsExclusivelyInGroup ${binPath}/halUpdateBranchLengths ${binPath}/halWriteNucleotides ${binPath}/halSetMetadata ${binPath}/halRenameGenomes ${binPath}/halRenameSequences
phastTargets = ${binPath}/ancestorsML ${binPath}/halAncestorsMLTest ${binPath}/ancestorsMLMP.py
ifdef ENABLE_PHYLOP
all : ${targets} ${phastTargets}
else
all : ${targets}
endif

clean : 
	rm -f ${targets} ${phastTargets} ./*.o

${binPath}/halRemoveGenome: halRemoveGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halRemoveGenome halRemoveGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halAddToBranch: halAddToBranch.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAddToBranch halAddToBranch.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halReplaceGenome: halReplaceGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halReplaceGenome halReplaceGenome.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/halAppendSubtree: halAppendSubtree.cpp markAncestors.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAppendSubtree halAppendSubtree.cpp markAncestors.o ${libPath}/halLib.a ${basicLibs}

${binPath}/findRegionsExclusivelyInGroup: findRegionsExclusivelyInGroup.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/findRegionsExclusivelyInGroup findRegionsExclusivelyInGroup.cpp ${libPath}/halLib.a ${libPath}/halLiftover.a ${basicLibs}

ancestorsMLBed.o: ancestorsMLBed.cpp
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -c $< -o $@ ${basicLibs} ${libPath}/halLib.a

ancestorsML.o: ancestorsML.cpp
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -c $< -o $@ ${basicLibs} ${libPath}/halLib.a

${binPath}/ancestorsML: ancestorsMLMain.cpp ancestorsML.o ancestorsMLBed.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -o ${binPath}/ancestorsML ancestorsMLMain.cpp ${basicLibs} ancestorsML.o ancestorsMLBed.o ${libPath}/halLib.a ${libPath}/halLiftover.a ${phastLinkflags}

${binPath}/ancestorsMLMP.py: ancestorsMLMP.py
	cp $< $@
	chmod +x $@

${binPath}/halAncestorsMLTest: ancestorsML.o ancestorsMLTest.cpp ancestorsMLBed.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} ${phastCflags} -o $@ ancestorsMLTest.cpp ${basicLibs} ancestorsML.o ${libPath}/halLib.a ${libPath}/halLiftover.a ${phastLinkflags}

${binPath}/halUpdateBranchLengths: halUpdateBranchLengths.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halUpdateBranchLengths halUpdateBranchLengths.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halWriteNucleotides: halWriteNucleotides.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halWriteNucleotides halWriteNucleotides.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halSetMetadata: halSetMetadata.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halSetMetadata halSetMetadata.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halRenameGenomes: halRenameGenomes.cpp renameFile.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o $@ $^

${binPath}/halRenameSequences: halRenameSequences.cpp renameFile.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o $@ $^

renameFile.o: renameFile.cpp
	${cpp} ${cppflags} -I ${libPath} -c -o $@ $<

markAncestors.o: markAncestors.cpp
	${cpp} ${cppflags} -I ${libPath} -c -o markAncestors.o markAncestors.cpp
