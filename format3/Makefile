rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/format3

# liftmap and twobit64, from sibling checkouts by default (see include.mk).  liftmap
# built with HTTP=1 needs libcurl.
inclSpec += -I${liftmapDir}/src -I${twobit64Dir}/src
format3Libs = ${libHal3} ${twobit64Dir}/libtwobit64.a ${liftmapDir}/libliftmap.a
format3LDLIBS = -lz -lpthread $(shell nm ${liftmapDir}/libliftmap.a 2>/dev/null | grep -q ' U curl_easy_init' && echo -lcurl)

# The core (manifest, validation) uses no HAL API and so no HDF5.  Only the converter
# links libHal.
libHal3_srcs = impl/hal3Manifest.cpp impl/hal3Alignment.cpp impl/hal3Liftover.cpp impl/hal3MafExport.cpp
libHal3_objs = ${libHal3_srcs:%.cpp=${modObjDir}/%.o}
halToFormat3_srcs = impl/halToFormat3Main.cpp impl/hal3Convert.cpp
halToFormat3_objs = ${halToFormat3_srcs:%.cpp=${modObjDir}/%.o}
hal3Validate_srcs = impl/hal3ValidateMain.cpp
hal3Validate_objs = ${hal3Validate_srcs:%.cpp=${modObjDir}/%.o}
hal3Liftover_srcs = impl/hal3LiftoverMain.cpp
hal3Liftover_objs = ${hal3Liftover_srcs:%.cpp=${modObjDir}/%.o}
hal3ToMaf_srcs = impl/hal3ToMafMain.cpp
hal3ToMaf_objs = ${hal3ToMaf_srcs:%.cpp=${modObjDir}/%.o}
# halLiftover's BED reading and writing, used as is (they need nothing of libHal but
# halCommon's string helpers), compiled here so no HAL library is linked
halShared_srcs = ${rootDir}/liftover/impl/halBedLine.cpp ${rootDir}/liftover/impl/halBedScanner.cpp ${rootDir}/api/impl/halCommon.cpp
halShared_objs = $(addprefix ${modObjDir}/shared/,$(notdir ${halShared_srcs:%.cpp=%.o}))
srcs = ${libHal3_srcs} ${halToFormat3_srcs} ${hal3Validate_srcs} ${hal3Liftover_srcs} ${hal3ToMaf_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halToFormat3 ${binDir}/hal3Validate ${binDir}/hal3Liftover ${binDir}/hal3ToMaf
otherLibs = ${format3Libs}
LDLIBS += ${format3LDLIBS}

# a plain C++ compiler for the core's programs: no h5c++ wrapper, so no HDF5
hal3CXX ?= c++

all : libs progs
libs: ${libHal3}
progs: ${progs}

${twobit64Dir}/libtwobit64.a ${liftmapDir}/libliftmap.a:
	@echo "format3 needs $@: build liftmap and twobit64 (sibling checkouts, or set liftmapDir / twobit64Dir)" >&2
	@false

${binDir}/hal3Validate: ${hal3Validate_objs} ${format3Libs}
	@mkdir -p $(dir $@)
	${hal3CXX} -std=c++11 -o $@ ${hal3Validate_objs} ${format3Libs} ${format3LDLIBS}

${binDir}/hal3Liftover: ${hal3Liftover_objs} ${halShared_objs} ${format3Libs}
	@mkdir -p $(dir $@)
	${hal3CXX} -std=c++11 -o $@ ${hal3Liftover_objs} ${halShared_objs} ${format3Libs} ${format3LDLIBS}

${binDir}/hal3ToMaf: ${hal3ToMaf_objs} ${halShared_objs} ${format3Libs}
	@mkdir -p $(dir $@)
	${hal3CXX} -std=c++11 -o $@ ${hal3ToMaf_objs} ${halShared_objs} ${format3Libs} ${format3LDLIBS}

${modObjDir}/shared/%.o: ${rootDir}/liftover/impl/%.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

${modObjDir}/shared/%.o: ${rootDir}/api/impl/%.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

clean :
	rm -f ${libHal3} ${objs} ${halShared_objs} ${progs} ${depends}

test: all
	./tests/format3Test.sh ${binDir}

include ${rootDir}/rules.mk

# header dependencies, written by rules.mk as each object is compiled
-include ${depends}
