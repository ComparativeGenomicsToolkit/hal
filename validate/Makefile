rootPath = ../
include ../include.mk

libSources = halValidateMain.cpp

all : ${binPath}/halValidate

clean : 
	rm -f ${binPath}/halValidate

${binPath}/halValidate : ${libSources} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halValidate ${libSources} ${libPath}/halLib.a ${basicLibs}

