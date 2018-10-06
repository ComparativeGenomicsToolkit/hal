rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/randgen

halRandGen_srcs = halRandGen.cpp
halRandGen_objs = ${halRandGen_srcs:%.cpp=${modObjDir}/%.o}
halTestGen_srcs = halTestGen.cpp
halTestGen_objs = ${halTestGen_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${halRandGen_srcs} ${halTestGen_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halRandGen ${binDir}/halTestGen

inclSpec += -I${halApiTestIncl}
otherLibs += ${halApiTestSupportLibs}

all: progs
libs:
progs: ${progs}

clean: 
	rm -f ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

