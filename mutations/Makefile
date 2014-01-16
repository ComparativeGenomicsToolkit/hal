rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources3=$(subst impl/halBranchMutationsMain.cpp,,${libSourcesAll})
libSources2=$(subst impl/halSummarizeMutationsMain.cpp,,${libSources3})
libSources1=$(subst impl/halSnps.cpp,,${libSources2})
libSources =$(subst impl/halIndels.cpp,,${libSources1})
libHeaders = $(wildcard inc/*.h)

all : ${libPath}/halMutations.a ${binPath}/halBranchMutations ${binPath}/halSummarizeMutations ${binPath}/halTreeMutations.py ${binPath}/halSnps ${binPath}/halIndels

clean : 
	rm -f ${libPath}/halMutations.a ${libPath}/*.h ${binPath}/halBranchMutations ${binPath}/halSummarizeMutations ${binPath}/halTreeMutations.py ${binPath}/halSnps ${binPath}/halIndels

${libPath}/halMutations.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMutations.a *.o
	ranlib halMutations.a 
	rm *.o
	mv halMutations.a ${libPath}/

${binPath}/halBranchMutations : impl/halBranchMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halBranchMutations impl/halBranchMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halSummarizeMutations : impl/halSummarizeMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halSummarizeMutations impl/halSummarizeMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halTreeMutations.py : impl/halTreeMutations.py
	cp impl/halTreeMutations.py ${binPath}
	chmod u+x ${binPath}/halTreeMutations.py

${binPath}/halSnps : impl/halSnps.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halSnps impl/halSnps.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halIndels : impl/halIndels.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halIndels impl/halIndels.cpp ${libPath}/halLib.a ${basicLibs}
