
rootPath = ../
include ../include.mk

libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${binPath}/halRemoveGenome ${binPath}/insertAlignedGenome ${binPath}/addNewLeafGenome ${binPath}/replaceGenome ${binPath}/halAppendSubtree

clean : 
	rm ${binPath}/halRemoveGenome ${binPath}/insertAlignedGenome ${binPath}/addNewLeafGenome

${binPath}/halRemoveGenome: halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halRemoveGenome halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/insertAlignedGenome: insertAlignedGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/insertAlignedGenome insertAlignedGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/addNewLeafGenome: addNewLeafGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/addNewLeafGenome addNewLeafGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/replaceGenome: replaceGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/replaceGenome replaceGenome.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halAppendSubtree: halAppendSubtree.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halAppendSubtree halAppendSubtree.cpp ${libPath}/halLib.a ${basicLibs}