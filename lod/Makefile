rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/lod

libHalLod_srcs = impl/halLodBlock.cpp impl/halLodExtract.cpp impl/halLodGraph.cpp \
    impl/halLodManager.cpp impl/halLodSegment.cpp
libHalLod_objs = ${libHalLod_srcs:%.cpp=${modObjDir}/%.o}
halLodExtract_srcs =impl/halLodExtractMain.cpp
halLodExtract_objs = ${halLodExtract_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalLod_srcs} ${halLodExtract_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
pyprogs = ${binDir}/halLodInterpolate.py
progs = ${binDir}/halLodExtract ${pyprogs}
otherLibs = ${libHalLod}

all : libs progs
libs: ${libHalLod}
progs: ${progs}

clean : 
	rm -f ${libHalLod} ${objs} ${progs} ${depends}
test:

${binDir}/%.py: %.py
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod +x,-w $@

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

