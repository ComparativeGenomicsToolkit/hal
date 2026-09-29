rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/format3

# liftmap and twobit64, from sibling checkouts by default (see include.mk).  liftmap
# built with HTTP=1 needs libcurl.
inclSpec += -I${liftmapDir}/src -I${twobit64Dir}/src
format3Libs = ${libHal3} ${twobit64Dir}/libtwobit64.a ${liftmapDir}/libliftmap.a
format3LDLIBS = -lz -lpthread $(shell nm ${liftmapDir}/libliftmap.a 2>/dev/null | grep -q ' U curl_easy_init' && echo -lcurl)

# One binary, bin/hal, with a command per file (impl/cmd_*.cpp).  Everything but
# `hal convert` uses no HAL API and no HDF5; convert reads the old format through
# libHal, so it is built in only with hal3Hdf5=1 (the default).  hal3Hdf5=0 builds a
# hal with no HDF5 at all.
hal3Hdf5 ?= 1

libHal3_srcs = impl/hal3Manifest.cpp impl/hal3Store.cpp impl/hal3Likelihood.cpp impl/hal3Alignment.cpp impl/hal3Liftover.cpp impl/hal3BedLiftover.cpp impl/hal3Walk.cpp \
    impl/hal3Universal.cpp impl/hal3MafExport.cpp impl/hal3Api.cpp impl/hal3Writer.cpp impl/hal3Cactus.cpp impl/hal3Update.cpp impl/hal3Cli.cpp
libHal3_objs = ${libHal3_srcs:%.cpp=${modObjDir}/%.o}
cmd_srcs = impl/cmd_validate.cpp impl/cmd_liftover.cpp impl/cmd_maf.cpp impl/cmd_fasta.cpp impl/cmd_stats.cpp impl/cmd_coverage.cpp impl/cmd_pack.cpp impl/cmd_update.cpp impl/cmd_merge.cpp impl/cmd_import.cpp impl/cmd_likelihoods.cpp
convert_srcs = impl/cmd_convert.cpp impl/hal3Convert.cpp
# halLiftover's BED reading and writing, used as is (they need nothing of libHal but
# halCommon's string helpers), compiled here so no HAL library is linked
halShared_srcs = ${rootDir}/liftover/impl/halBedLine.cpp ${rootDir}/liftover/impl/halBedScanner.cpp ${rootDir}/api/impl/halCommon.cpp
halShared_objs = $(addprefix ${modObjDir}/shared/,$(notdir ${halShared_srcs:%.cpp=%.o}))

ifeq (${hal3Hdf5},1)
hal_objs = ${modObjDir}/impl/hal_convert.o ${cmd_srcs:%.cpp=${modObjDir}/%.o} ${convert_srcs:%.cpp=${modObjDir}/%.o}
else
hal_objs = ${modObjDir}/impl/hal_core.o ${cmd_srcs:%.cpp=${modObjDir}/%.o}
endif

srcs = ${libHal3_srcs} ${cmd_srcs} ${convert_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o} ${modObjDir}/impl/hal_convert.o ${modObjDir}/impl/hal_core.o
depends = ${srcs:%.cpp=%.depend} impl/hal.depend
progs = ${binDir}/hal
otherLibs = ${format3Libs}
LDLIBS += ${format3LDLIBS}

# a plain C++ compiler for linking without HDF5 (no h5c++ wrapper)
hal3CXX ?= c++

all : libs progs
libs: ${libHal3}
progs: ${progs}

${twobit64Dir}/libtwobit64.a ${liftmapDir}/libliftmap.a:
	@echo "format3 needs $@: build liftmap and twobit64 (sibling checkouts, or set liftmapDir / twobit64Dir)" >&2
	@false

# the command table, with and without convert
${modObjDir}/impl/hal_convert.o: impl/hal.cpp
	@mkdir -p $(dir $@)
	${CXX} -MM -MT $@ ${CXXFLAGS} ${inclSpec} -c $< >impl/hal.depend
	${CXX} ${CXXFLAGS} ${inclSpec} -DHAL3_CONVERT -c $< -o $@
${modObjDir}/impl/hal_core.o: impl/hal.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

ifeq (${hal3Hdf5},1)
${binDir}/hal: ${hal_objs} ${halShared_objs} ${format3Libs} ${libHal}
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} -o $@ ${hal_objs} ${halShared_objs} ${format3Libs} ${libHal} ${LDLIBS}
else
${binDir}/hal: ${hal_objs} ${halShared_objs} ${format3Libs}
	@mkdir -p $(dir $@)
	${hal3CXX} -std=c++11 -o $@ ${hal_objs} ${halShared_objs} ${format3Libs} ${format3LDLIBS}
endif

${modObjDir}/shared/%.o: ${rootDir}/liftover/impl/%.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

${modObjDir}/shared/%.o: ${rootDir}/api/impl/%.cpp
	@mkdir -p $(dir $@)
	${CXX} ${CXXFLAGS} ${inclSpec} -c $< -o $@

clean :
	rm -f ${libHal3} ${objs} ${halShared_objs} ${progs} ${depends} ${binDir}/hal3PackTest

${binDir}/hal3PackTest: tests/packTest.cpp ${format3Libs}
	@mkdir -p $(dir $@)
	${hal3CXX} -std=c++11 ${inclSpec} -o $@ $< ${format3Libs} ${format3LDLIBS}

test: all ${binDir}/hal3PackTest
	./tests/format3Test.sh ${binDir}

include ${rootDir}/rules.mk

# header dependencies, written by rules.mk as each object is compiled
-include ${depends}
