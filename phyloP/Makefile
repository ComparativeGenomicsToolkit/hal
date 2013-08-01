rootPath = ../
include ../include.mk

libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/halPhyloPMain.cpp,,${libSourcesAll})
libHeaders = $(wildcard inc/*.h)

ifdef ENABLE_PHYLOP
all : ${libPath}/halPhyloP.a ${binPath}/halPhyloP ${binPath}/halPhyloPTrain.py ${binPath}/halPhyloPMP.py ${binPath}/halTreePhyloP.py
else
all :
	cd .
endif

clean : 
	rm -f ${libPath}/halPhyloP.a ${libPath}/halPhyloP.h ${binPath}/halPhyloP *.o ${binPath}/halPhyloPTrain.py ${binPath}/halPhyloPMP.py ${binPath}/halTreePhyloP.py

${libPath}/halPhyloP.a : ${libSources} ${libHeaders} ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	rm -f *.o
	${cpp} ${phyloPcppflags} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halPhyloP.a *.o
	ranlib halPhyloP.a 
	rm *.o
	mv halPhyloP.a ${libPath}/

${binPath}/halPhyloP : impl/halPhyloPMain.cpp ${libPath}/halPhyloP.a ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${phyloPcppflags} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halPhyloP impl/halPhyloPMain.cpp ${libPath}/halPhyloP.a ${libPath}/halLib.a ${libPath}/halLiftover.a ${basicLibs} ${phyloPlibs}

${binPath}/halPhyloPTrain.py : halPhyloPTrain.py
	cp halPhyloPTrain.py ${binPath}/halPhyloPTrain.py
	chmod +x ${binPath}/halPhyloPTrain.py

${binPath}/halPhyloPMP.py : halPhyloPMP.py
	cp halPhyloPMP.py ${binPath}/halPhyloPMP.py
	chmod +x ${binPath}/halPhyloPMP.py

${binPath}/halTreePhyloP.py : halTreePhyloP.py
	cp halTreePhyloP.py ${binPath}/halTreePhyloP.py
	chmod +x ${binPath}/halTreePhyloP.py
