rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/halPhyloPMain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)

all : ${libPath}/halPhyloP.a ${binPath}/halPhyloP ${binPath}/halPhyloPTrain.py

clean : 
	rm -f ${libPath}/halPhyloP.a ${libPath}/halPhyloP.h ${binPath}/halPhyloP *.o ${binPath}/halPhyloPTrain.py

${libPath}/halPhyloP.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources} ${phyloPcppflags}
	ar rc halPhyloP.a *.o
	ranlib halPhyloP.a 
	rm *.o
	mv halPhyloP.a ${libPath}/

${binPath}/halPhyloP : impl/halPhyloPMain.cpp ${libPath}/halPhyloP.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halPhyloP impl/halPhyloPMain.cpp ${libPath}/halPhyloP.a ${libPath}/halLib.a ${basicLibs} ${phyloPcppflags} ${phyloPlibs}

${binPath}/halPhyloPTrain.py : halPhyloPTrain.py
	cp halPhyloPTrain.py ${binPath}/halPhyloPTrain.py
	chmod +x ${binPath}/halPhyloPTrain.py
