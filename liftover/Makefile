rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/halLiftoverMain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halLiftover.a ${binPath}/halLiftover

clean : 
	rm -f ${libPath}/halLiftover.a ${libPath}/*.h ${binPath}/halLiftover 

${libPath}/halLiftover.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halLiftover.a *.o
	ranlib halLiftover.a 
	rm *.o
	mv halLiftover.a ${libPath}/

${binPath}/halLiftover : impl/halLiftoverMain.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halLiftover impl/halLiftoverMain.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

