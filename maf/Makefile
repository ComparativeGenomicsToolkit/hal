rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources1=$(subst impl/hal2maf.cpp,,${libSourcesAll})
libSources=$(subst impl/maf2hal.cpp,,${libSources1})
libHeaders = $(wildcard inc/*.h )
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/allTests.h
libTests = $(wildcard tests/*.cpp)
libTestsHeaders = $(wildcard tests/*.h)
libHalTestsAll := $(wildcard ../api/tests/*.cpp)
libHalTests = $(subst ../api/tests/allTests.cpp,,${libHalTestsAll})

all : ${libPath}/halMaf.a ${binPath}/hal2maf ${binPath}/maf2hal ${binPath}/halMafTests ${binPath}/hal2mafMP.py

clean : 
	rm -f ${libPath}/halMaf.a ${libPath}/*.h ${binPath}/hal2maf ${binPath}/maf2hal ${binPath}/hal2mafMP.py

${libPath}/halMaf.a : ${libSources} ${libHeaders} ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMaf.a *.o
	ranlib halMaf.a 
	rm *.o
	mv halMaf.a ${libPath}/

${binPath}/hal2maf : impl/hal2maf.cpp ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/hal2maf impl/hal2maf.cpp ${libPath}/halMaf.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

${binPath}/maf2hal : impl/maf2hal.cpp ${libPath}/halMaf.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/maf2hal impl/maf2hal.cpp ${libPath}/halMaf.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halMafTests : ${libTests} ${libTestsHeaders} ${libTestsCommon} ${libTestsHeadersCommon} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -I ../api/tests -o ${binPath}/halMafTests ${libHalTests} ${libTests} ${libPath}/halLib.a ${libPath}/halMaf.a ${basicLibs}

${binPath}/hal2mafMP.py : hal2mafMP.py
	cp hal2mafMP.py ${binPath}/hal2mafMP.py
	chmod +x ${binPath}/hal2mafMP.py
