rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/maf

libHalMaf_srcs = impl/halMafBed.cpp impl/halMafBlock.cpp impl/halMafExport.cpp \
    impl/halMafScanDimensions.cpp impl/halMafScanner.cpp impl/halMafScanReference.cpp \
    impl/halMafWriteGenomes.cpp
libHalMaf_objs = ${libHalMaf_srcs:%.cpp=${modObjDir}/%.o}
hal2maf_srcs = impl/hal2maf.cpp
hal2maf_objs = ${hal2maf_srcs:%.cpp=${modObjDir}/%.o}
maf2hal_srcs = impl/maf2hal.cpp
maf2hal_objs = ${maf2hal_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalMaf_srcs} ${hal2maf_srcs} ${maf2hal_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/hal2maf ${binDir}/maf2hal
otherLibs = ${libHalMaf} ${libHalLiftover}
inclSpec += -I${rootDir}/liftover/inc

all: libs progs
libs: ${libHalMaf}
progs: ${progs}

clean : 
	rm -f ${libHalMaf} ${objs} ${progs} ${depends}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

