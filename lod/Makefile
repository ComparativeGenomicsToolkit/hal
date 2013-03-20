rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/halLodExtractMain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halLod.a ${binPath}/halLodExtract

clean : 
	rm -f ${libPath}/halLod.a ${libPath}/*.h ${binPath}/halLodExtract 

${libPath}/halLod.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halLod.a *.o
	ranlib halLod.a 
	rm *.o
	mv halLod.a ${libPath}/

${binPath}/halLodExtract : impl/halLodExtractMain.cpp ${libPath}/halLod.a ${libPath}/halLod.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halLodExtract impl/halLodExtractMain.cpp ${libPath}/halLod.a ${libPath}/halLib.a ${basicLibs}
