rootPath = ../
include ../include.mk

libSources = hal2fasta.cpp 

all : ${binPath}/hal2fasta

clean : 
	rm -f ${binPath}/hal2fasta

${binPath}/hal2fasta : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/hal2fasta ${libSources} ${libPath}/halLib.a ${basicLibs}

