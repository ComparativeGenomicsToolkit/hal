rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/halMaskMain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halMask.a ${binPath}/halMask

clean : 
	rm -f ${libPath}/halMask.a ${libPath}/*.h ${binPath}/halMask 

${libPath}/halMask.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMask.a *.o
	ranlib halMask.a 
	rm *.o
	mv halMask.a ${libPath}/

${binPath}/halMask : impl/halMaskMain.cpp ${libPath}/halMask.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halMask impl/halMaskMain.cpp ${libPath}/halMask.a ${libPath}/halLib.a ${basicLibs}

