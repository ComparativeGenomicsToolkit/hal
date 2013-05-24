rootPath = ../
include ../include.mk

libRandSources = halRandGen.cpp ${rootPath}/api/tests/halRandomData.cpp
libTestGenSources = halTestGen.cpp ${rootPath}/api/tests/halAlignmentTest.cpp ${rootPath}/api/tests/halAlignmentInstanceTest.cpp ${rootPath}/api/tests/halRearrangementTest.cpp ${rootPath}/api/tests/halTopSegmentTest.cpp ${rootPath}/api/tests/halBottomSegmentTest.cpp

all : ${binPath}/halRandGen ${binPath}/halTestGen

clean : 
	rm -f ${binPath}/halRandGen ${binPath}/halTestGen

${binPath}/halRandGen : Makefile ${libRandSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halRandGen ${libRandSources} ${libPath}/halLib.a ${basicLibs}

${binPath}/halTestGen : Makefile ${libTestGenSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halTestGen ${libTestGenSources} ${libPath}/halLib.a ${basicLibs}