rootPath = ../
include ../include.mk

libSources = halRandGen.cpp ${rootPath}/api/tests/halRandomData.cpp

all : ${binPath}/halRandGen

clean : 
	rm -f ${binPath}/halRandGen

${binPath}/halRandGen : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halRandGen ${libSources} ${libPath}/halLib.a ${basicLibs}

