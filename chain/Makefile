rootPath = ../
include ../include.mk

libSourcesAll = impl/*.cpp
libSources=$(subst impl/hal2chain.cpp,,${libSourcesAll})
libHeaders = inc/*.h 

all : ${libPath}/halChain.a ${binPath}/hal2chain test/blockVizTest test/blockInterpolateTest

clean : 
	rm -f ${libPath}/halChain.a ${libPath}/*.h ${binPath}/hal2chain test/blockVizTest test/blockInterpolateTest

${libPath}/halChain.a : ${libSources} ${libHeaders} ${libPath}/halLib.a ${basicLibsDependencies} 
	cp ${libHeaders} ${libPath}/
	${cpp} ${cppflags} -I inc -I hdf5_impl -I impl -I ${libPath}/ -c ${libSources}
	ar rc halChain.a *.o
	ranlib halChain.a 
	rm *.o
	mv halChain.a ${libPath}/

${binPath}/hal2chain : impl/hal2chain.cpp ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/hal2chain impl/hal2chain.cpp ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibs}

test/blockVizTest: test/blockVizTest.c ${libPath}/halChain.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockVizTest test/blockVizTest.c ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibs}
 
test/blockInterpolateTest: test/blockInterpolateTest.cpp ${libPath}/halChain.a ${libPath}/halLib.a ${basicLib\
sDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o test/blockInterpolateTest test/blockInterpolateTest.cpp ${libPath}/halChain.a ${libPath}/halLib.a ${basicLibs}
