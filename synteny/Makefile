rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/synteny

hal2synteny_srcs = impl/main.cpp impl/hal2psl.cpp impl/psl_io.cpp impl/psl_merger.cpp
hal2synteny_objs = ${hal2synteny_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${hal2synteny_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
inclSpec += -I${rootDir}/liftover/inc
otherLibs += ${libHalLiftover}
progs = ${binDir}/hal2synteny

all: progs
libs:
progs: ${progs}

clean : 
	rm -f ${objs} ${progs} ${depends}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

