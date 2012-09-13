rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/hal2maf.cpp,,${libSourcesAll})
libHeaders = inc/*.h 
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = tests/*.cpp 
libTestsHeaders = tests/*.h 
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${libPath}/halMafExport.a ${binPath}/hal2maf ${binPath}/halMafTests

clean : 
	rm -f ${libPath}/halMafExport.a ${libPath}/*.h ${binPath}/hal2maf 

${libPath}/halMafExport.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMafExport.a *.o
	ranlib halMafExport.a 
	rm *.o
	mv halMafExport.a ${libPath}/

${binPath}/hal2maf : impl/hal2maf.cpp ${libPath}/halMafExport.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/hal2maf impl/hal2maf.cpp ${libPath}/halMafExport.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halMafTests : ${libTests} ${libTestsHeaders} ${libTestsCommon} ${libTestsHeadersCommon} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -I ../api/tests -o ${binPath}/halMafTests ${libHalTests} ${libTests} ${libPath}/halLib.a ${libPath}/halMafExport.a ${basicLibs}
