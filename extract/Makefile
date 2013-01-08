rootPath = ../
include ../include.mk

libSources = halExtract.cpp 
libAlignedSources = halAlignedExtract.cpp

all : ${binPath}/halExtract ${binPath}/halAlignedExtract

clean : 
	rm -f ${binPath}/halExtract ${binPath}/halAlignedExtract

${binPath}/halExtract : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halExtract ${libSources} ${libPath}/halLib.a ${basicLibs}

${binPath}/halAlignedExtract : ${libAlignedSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halAlignedExtract ${libAlignedSources} ${libPath}/halLib.a ${basicLibs}
