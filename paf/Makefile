rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/paf

hal2paf_srcs = hal2paf.cpp 
hal2paf_objs = ${hal2paf_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${hal2paf_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/hal2paf

all: progs
libs:
progs: ${progs}

clean: 
	rm -f  ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

