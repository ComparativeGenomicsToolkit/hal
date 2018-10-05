rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/stats

libHalStats_srcs = impl/halStats.cpp
libHalStats_objs = ${libHalStats_srcs:%.cpp=${modObjDir}/%.o}
halStats_srcs = impl/halStatsMain.cpp
halStats_objs = ${halStats_srcs:%.cpp=${modObjDir}/%.o}
halCoverage_srcs = impl/halCoverage.cpp
halCoverage_objs = ${halCoverage_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalStats_srcs} ${halStats_srcs} ${halCoverage_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halStats ${binDir}/halCoverage
otherLibs = ${libHalStats}

all : libs progs
libs: ${libHalStats}
progs: ${progs}

clean : 
	rm -f ${libHalStats} ${objs} ${progs} ${depends}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

