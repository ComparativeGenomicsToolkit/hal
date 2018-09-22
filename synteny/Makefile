rootPath = ../
include ../include.mk

libHeaders = $(wildcard inc/*)
libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/main.cpp,,${libSourcesAll})

cpp = h5c++ ${h5prefix}

all : ${binPath}/hal2synteny    

${binPath}/hal2synteny : impl/main.cpp  ${libPath}/halLib.a  ${libPath}/halLiftover.a  ${basicLibsDependencies}
	${cpp} ${cppflags} -std=c++11 -I${sonLibPath} -I${rootPath}/liftover/inc -I${rootPath}/api/inc -Iinc -Iimpl impl/main.cpp ${libPath}/halLib.a ${libPath}/halLiftover.a  ${sonLibPath}/sonLib.a -o ${binPath}/hal2synteny

clean:
	rm -f ${binPath}/hal2synteny
