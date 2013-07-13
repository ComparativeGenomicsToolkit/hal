rootPath = ../
include ../include.mk

libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${binPath}/halRemoveGenome ${binPath}/insertAlignedGenome ${binPath}/addNewLeafGenome

clean : 
	rm ${binPath}/halRemoveGenome ${binPath}/insertAlignedGenome ${binPath}/addNewLeafGenome copyGenomes.o

${binPath}/halRemoveGenome: halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/halRemoveGenome halRemoveGenome.cpp ${libPath}/halLib.a ${basicLibs}

copyGenomes.o: copyGenomes.cpp
	${cpp} -c ${cppflags} -I ${libPath} -o copyGenomes.o copyGenomes.cpp

${binPath}/insertAlignedGenome: insertAlignedGenome.cpp copyGenomes.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/insertAlignedGenome insertAlignedGenome.cpp ${libPath}/halLib.a ${basicLibs} copyGenomes.o

${binPath}/addNewLeafGenome: addNewLeafGenome.cpp copyGenomes.o ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I ${libPath} -o ${binPath}/addNewLeafGenome addNewLeafGenome.cpp ${libPath}/halLib.a ${basicLibs} copyGenomes.o
