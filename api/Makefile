rootPath = ../
include ../include.mk

libSources = impl/*.cpp hdf5_impl/*.cpp
libHeaders = inc/*.h 
libInternalHeaders = hdf5_impl/*.h impl/*.h
libTests = tests/*.cpp
libTestsHeaders = tests/*.h 
libHdf5Tests = hdf5_tests/*.cpp
libHdf5TestsHeaders = hdf5_tests/*.h

all : ${libPath}/halLib.a ${binPath}/halApiTests ${binPath}/halHdf5Tests 

clean : 
	rm -f ${libPath}/halLib.a ${libPath}/*.h ${binPath}/halApiTests ${binPath}/halHdf5Tests

doxy :
	doxygen doc/doxy.cfg

${libPath}/halLib.a : ${libSources} ${libHeaders} ${libInternalHeaders} ${basicLibsDependencies}
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halLib.a *.o
	ranlib halLib.a 
	rm *.o
	mv halLib.a ${libPath}/

${binPath}/halApiTests : ${libTests} ${libTestsHeaders} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halApiTests ${libTests} ${libPath}/halLib.a ${basicLibs}

${binPath}/halHdf5Tests : ${libHdf5Tests} ${libHdf5TestsHeaders} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I hdf5_impl -I ${libPath} -I impl -I tests -o ${binPath}/halHdf5Tests ${libHdf5Tests} ${libPath}/halLib.a ${basicLibs}
