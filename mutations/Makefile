rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/halMutationsMain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halMutations.a ${binPath}/halMutations

clean : 
	rm -f ${libPath}/halMutations.a ${libPath}/*.h ${binPath}/halMutations 

${libPath}/halMutations.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMutations.a *.o
	ranlib halMutations.a 
	rm *.o
	mv halMutations.a ${libPath}/

${binPath}/halMutations : impl/halMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halMutations impl/halMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibs}

