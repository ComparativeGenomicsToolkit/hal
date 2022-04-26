rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/blockViz

libHalBlockViz_srcs = impl/halBlockViz.cpp
libHalBlockViz_objs = ${libHalBlockViz_srcs:%.cpp=${modObjDir}/%.o}
blockVizBed_srcs = tests/blockVizBed.cpp
blockVizBed_objs = ${blockVizBed_srcs:%.cpp=${modObjDir}/%.o}
blockVizMaf_srcs = tests/blockVizMaf.cpp
blockVizMaf_objs = ${blockVizMaf_srcs:%.cpp=${modObjDir}/%.o}
blockVizTest_srcs = tests/blockVizTest.cpp
blockVizTest_objs = ${blockVizTest_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalBlockViz_srcs} ${blockVizBed_srcs} \
    ${blockVizMaf_srcs} ${blockVizTest_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
inclSpec += -I${rootDir}/liftover/inc -I${rootDir}/lod/inc -I${rootDir}/maf/inc -I${halApiTestIncl}
otherLibs += ${halApiTestSupportLibs} ${libHalBlockViz} ${libHalLiftover} ${libHalLod} ${libHalMaf}
progs =  ${binDir}/blockVizBed ${binDir}/blockVizMaf ${binDir}/blockVizTest

testTmpDir = output
testHdf5Hal = ${testTmpDir}/small.haf5.hal
testMmapHal = ${testTmpDir}/small.mmap.hal

all: libs progs
libs: ${libHalBlockViz}
progs: ${progs}

# only blockVizTest/Maf.cpp should include kent .h files, to preven conflicts with older versions of
# phast
${modObjDir}/tests/blockVizTest.o : tests/blockVizTest.cpp
	@mkdir -p $(dir $@)
	${CXX} -MM -MT $@ ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< >$*.depend
	${CXX} ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< -o $@

${modObjDir}/tests/blockVizMaf.o : tests/blockVizMaf.cpp
	@mkdir -p $(dir $@)
	${CXX} -MM -MT $@ ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< >$*.depend
	${CXX} ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< -o $@

${modObjDir}/tests/blockVizBed.o : tests/blockVizBed.cpp
	@mkdir -p $(dir $@)
	${CXX} -MM -MT $@ ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< >$*.depend
	${CXX} ${CXXFLAGS} ${UDCCXXFLAGS} ${inclSpec} -c $< -o $@

clean: 
	rm -f ${libHalBlockViz} ${objs} ${progs} ${depends}
	rm -rf ${testTmpDir}

test: blockVizHdf5Tests blockVizMmapTests

blockVizHdf5Tests: ${testHdf5Hal} ${progs}
	${binDir}/blockVizTest --verbose --doSeq ${testHdf5Hal} Genome_2 Genome_0 Genome_0_seq 0 3000 >${testTmpDir}/$@.out
	diff tests/expected/$@.out ${testTmpDir}/$@.out


blockVizMmapTests: ${testMmapHal} ${progs}
	${binDir}/blockVizTest --verbose --doSeq ${testMmapHal} Genome_2 Genome_0 Genome_0_seq 0 3000 >${testTmpDir}/$@.out
	diff tests/expected/$@.out ${testTmpDir}/$@.out

randGenArgs = --preset small --seed 0 --minSegmentLength 3000  --maxSegmentLength 5000

${testHdf5Hal}: ${progs} ${binDir}/halRandGen
	@mkdir -p $(dir $@)
	${binDir}/halRandGen ${randGenArgs} --format hdf5 $@

${testMmapHal}: ${progs} ${binDir}/halRandGen
	@mkdir -p $(dir $@)
	${binDir}/halRandGen ${randGenArgs} --format mmap $@

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

