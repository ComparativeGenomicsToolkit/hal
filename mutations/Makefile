rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources1=$(subst impl/halBranchMutationsMain.cpp,,${libSourcesAll})
libSources=$(subst impl/halSummarizeMutationsMain.cpp,,${libSources1})
libHeaders = inc/*.h 

all : ${libPath}/halMutations.a ${binPath}/halBranchMutations ${binPath}/halSummarizeMutations

clean : 
	rm -f ${libPath}/halMutations.a ${libPath}/*.h ${binPath}/halBranchMutations ${binPath}/halSummarizeMutations

${libPath}/halMutations.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halMutations.a *.o
	ranlib halMutations.a 
	rm *.o
	mv halMutations.a ${libPath}/

${binPath}/halBranchMutations : impl/halBranchMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halBranchMutations impl/halBranchMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halSummarizeMutations : impl/halSummarizeMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halSummarizeMutations impl/halSummarizeMutationsMain.cpp ${libPath}/halMutations.a ${libPath}/halLib.a ${basicLibs}
