rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/halStatsMain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)

all : ${libPath}/halStats.a ${binPath}/halStats

clean : 
	rm -f ${libPath}/halStats.a ${libPath}/*.h ${binPath}/halStats 

${libPath}/halStats.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halStats.a *.o
	ranlib halStats.a 
	rm *.o
	mv halStats.a ${libPath}/

${binPath}/halStats : impl/halStatsMain.cpp ${libPath}/halStats.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halStats impl/halStatsMain.cpp ${libPath}/halStats.a ${libPath}/halLib.a ${basicLibs}
	rm *.o
