rootPath = ../
include ../include.mk

libSources = impl/*.cpp hdf5_impl/*.cpp
libHeaders = inc/*.h 
libInternalHeaders = hdf5_impl/*.h
libTests = hdf5_tests/*.cpp
libTestsHeaders = hdf5_tests/*.h

all : ${libPath}/halLib.a ${binPath}/halTests 

clean : 
	rm -f ${libPath}/halLib.a ${libPath}/*.h ${binPath}/halTests

doxy :
	doxygen doc/doxy.cfg

${libPath}/halLib.a : ${libSources} ${libHeaders} ${libInternalHeaders} ${basicLibsDependencies}
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I ${libPath}/ -c ${libSources}
	ar rc halLib.a *.o
	ranlib halLib.a 
	rm *.o
	mv halLib.a ${libPath}/

${binPath}/halTests : ${libTests} ${libTestsHeaders} ${libSources} ${libHeaders} ${libInternalHeaders} ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I hdf5_impl -I ${libPath} -I impl -I tests -o ${binPath}/halTests ${libTests} ${libPath}/halLib.a ${basicLibs}
