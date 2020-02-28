rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/synteny

halSynteny_srcs = impl/halSynteny.cpp impl/hal2psl.cpp impl/psl_io.cpp impl/psl_merger.cpp
halSynteny_objs = ${halSynteny_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${halSynteny_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
inclSpec += -I${rootDir}/liftover/inc
otherLibs += ${libHalLiftover}
progs = ${binDir}/halSynteny

all: progs
libs:
progs: ${progs}

clean : 
	rm -f ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:
