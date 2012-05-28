rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/halConsMain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halCons.a ${binPath}/halCons

clean : 
	rm -f ${libPath}/halCons.a ${libPath}/*.h ${binPath}/halCons 

${libPath}/halCons.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halCons.a *.o
	ranlib halCons.a 
	rm *.o
	mv halCons.a ${libPath}/

${binPath}/halCons : impl/halConsMain.cpp ${libPath}/halCons.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halCons impl/halConsMain.cpp ${libPath}/halCons.a ${libPath}/halLib.a ${basicLibs}

