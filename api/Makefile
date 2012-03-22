rootPath = ../
include ../include.mk

libSources = impl/*.cpp
libHeaders = inc/*.h
libInternalHeaders = impl/*.h
libTests = tests/*.cpp
libTestsHeaders = tests/*.h

all : ${libPath}/halLib.a ${binPath}/halTests 

clean : 
	rm -f ${libPath}/halLib.a ${libPath}/*.h ${binPath}/halTests

${libPath}/halLib.a : ${libSources} ${libHeaders} ${libInternalHeaders} ${basicLibsDependencies}
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc halLib.a *.o
	ranlib halLib.a 
	rm *.o
	mv halLib.a ${libPath}/

${binPath}/halTests : ${libTests} ${libTestsHeaders} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I ${libPath} -I impl -I tests -o ${binPath}/halTests ${libTests} ${libPath}/halLib.a ${basicLibs}
