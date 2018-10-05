rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/alignmentDepth

halAlignmentDepth_srcs = halAlignmentDepth.cpp
halAlignmentDepth_objs = ${halAlignmentDepth_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${halAlignmentDepth_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halAlignmentDepth

all: progs
libs:
progs: ${progs}

clean: 
	rm -f ${objs} ${progs} ${depends}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

