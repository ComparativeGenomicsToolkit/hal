rootPath = ../
include ../include.mk

libSources = halAlignmentDepth.cpp 

all : ${binPath}/halAlignmentDepth

clean : 
	rm -f ${binPath}/halAlignmentDepth

${binPath}/halAlignmentDepth : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	rm -f ${binPath}/halAlignability
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halAlignmentDepth ${libSources} ${libPath}/halLib.a ${basicLibs}

