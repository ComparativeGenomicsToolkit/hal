rootPath = ../
include ../include.mk

libSources = halAlignability.cpp 

all : ${binPath}/halAlignability

clean : 
	rm -f ${binPath}/halAlignability

${binPath}/halAlignability : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halAlignability ${libSources} ${libPath}/halLib.a ${basicLibs}

