rootPath = ../
include ../include.mk

libHeaders = $(wildcard inc/*)
libSourcesAll = $(wildcard impl/*.cpp)
libSources=$(subst impl/main.cpp,,${libSourcesAll})

cpp = h5c++ ${h5prefix}

all : ${binPath}/halSynteny    

${binPath}/halSynteny : impl/main.cpp  ${libPath}/halLib.a  ${libPath}/halLiftover.a  ${basicLibsDependencies}
	${cpp} ${cppflags} -std=c++11 -I${sonLibPath} -I${rootPath}/liftover/inc -I${rootPath}/api/inc -Iinc -Iimpl impl/main.cpp ${libPath}/halLib.a ${libPath}/halLiftover.a  ${sonLibPath}/sonLib.a -o ${binPath}/halSynteny

clean:
	rm -f ${binPath}/halSynteny
