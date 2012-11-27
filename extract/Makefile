rootPath = ../
include ../include.mk

libSources = halExtract.cpp 

all : ${binPath}/halExtract

clean : 
	rm -f ${binPath}/halExtract

${binPath}/halExtract : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halExtract ${libSources} ${libPath}/halLib.a ${basicLibs}

