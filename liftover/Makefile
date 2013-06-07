rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/halLiftoverMain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)
libTestSources = $(wildcard tests/*.cpp)
libTestHeaders = $(wildcard tests/*.h)
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h

all : ${libPath}/halLiftover.a ${binPath}/halLiftover ${binPath}/halLiftoverTests

clean : 
	rm -f ${libPath}/halLiftover.a ${libPath}/*.h ${binPath}/halLiftover ${binPath}/halLiftoverTests

${libPath}/halLiftover.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halLiftover.a *.o
	ranlib halLiftover.a 
	rm *.o
	mv halLiftover.a ${libPath}/

${binPath}/halLiftover : impl/halLiftoverMain.cpp ${libPath}/halLiftover.a ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halLiftover impl/halLiftoverMain.cpp ${libPath}/halLiftover.a ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halLiftoverTests : ${libTestSources} ${libTestHeaders} ${libTestsCommon} ${libTestsHeadersCommon} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -I ../api/tests -o ${binPath}/halLiftoverTests  ${libTestSources} ${libTestsCommon} ${libPath}/halChain.a ${libPath}/halLib.a ${libPath}/halLiftover.a ${basicLibs}