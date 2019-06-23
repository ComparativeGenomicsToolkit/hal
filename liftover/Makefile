rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/liftover

libHalLiftover_srcs = impl/halBedLine.cpp impl/halBedScanner.cpp impl/halBlockLiftover.cpp \
    impl/halBlockMapper.cpp impl/halColumnLiftover.cpp impl/halLiftover.cpp \
    impl/halWiggleLiftover.cpp impl/halWiggleLoader.cpp impl/halWiggleScanner.cpp
libHalLiftover_objs = ${libHalLiftover_srcs:%.cpp=${modObjDir}/%.o}
halLiftover_srcs = impl/halLiftoverMain.cpp
halLiftover_objs = ${halLiftover_srcs:%.cpp=${modObjDir}/%.o}
halWiggleLiftover_srcs = impl/halWiggleLiftoverMain.cpp
halWiggleLiftover_objs = ${halWiggleLiftover_srcs:%.cpp=${modObjDir}/%.o}
halLiftoverTests_srcs = tests/halLiftoverTests.cpp
halLiftoverTests_objs = ${halLiftoverTests_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalLiftover_srcs} ${halLiftover_srcs} ${halWiggleLiftover_srcs} ${halLiftover_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halLiftover ${binDir}/halWiggleLiftover ${binDir}/halLiftoverTests
otherLibs = ${libHalLiftover}

# tests use api/tests/halAlignmentTest
inclSpec += -I${rootDir}/api/tests
otherLibs += ${rootDir}/objs/api/tests/halAlignmentTest.o ${rootDir}/objs/api/tests/halAlignmentInstanceTest.o

all : libs progs
libs: ${libHalLiftover}
progs: ${progs}

clean : 
	rm -f ${libHalLiftover} ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

