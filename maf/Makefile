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
halMafTests_srcs = tests/halMafTests.cpp tests/halMafBlockTest.cpp tests/halMafExportTest.cpp
halMafTests_objs = ${halMafTests_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalMaf_srcs} ${hal2maf_srcs} ${maf2hal_srcs} ${halMafTests_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
pyprogs = ${binDir}/hal2mafMP.py
progs = ${binDir}/hal2maf ${binDir}/maf2hal ${binDir}/halMafTests ${pyprogs}
otherLibs = ${libHalMaf} ${libHalLiftover} ${halApiTestSupportLibs}
inclSpec += -I${rootDir}/liftover/inc -I${halApiTestIncl}

all: libs progs
libs: ${libHalMaf}
progs: ${progs}

clean : 
	rm -rf ${libHalMaf} ${objs} ${progs} ${depends} output

test: halMafTests hal2mafCmdTests naiveLiftUpTests

halMafTests:
	${binDir}/halMafTests

naiveLiftUpTests:
	python2 -m pytest impl/naiveLiftUp.py

hal2mafCmdTests: hal2MafSmallMMapTest hal2MafSmallHdf5Test hal2MafSeqTest hal2MafSeqPartTest

hal2MafSmallMMapTest: output/small.mmap.hal
	../bin/hal2maf output/small.mmap.hal output/$@.maf
	diff tests/expected/hal2MafSmallTest.maf output/$@.maf

hal2MafSmallHdf5Test: output/small.hdf5.hal
	../bin/hal2maf output/small.hdf5.hal output/$@.maf
	diff tests/expected/hal2MafSmallTest.maf output/$@.maf

hal2MafSeqTest: output/small.mmap.hal
	../bin/hal2maf --refGenome Genome_2 --refSequence Genome_2_seq output/small.mmap.hal output/$@.maf
	diff tests/expected/$@.maf output/$@.maf

hal2MafSeqPartTest: output/small.mmap.hal
	../bin/hal2maf --refGenome Genome_2 --refSequence Genome_2_seq --start 1000 --length 2000 output/small.mmap.hal output/$@.maf
	diff tests/expected/$@.maf output/$@.maf

output/small.mmap.hal:
	@mkdir -p output
	../bin/halRandGen --preset small --seed 0 --testRand --format mmap output/small.mmap.hal

output/small.hdf5.hal:
	@mkdir -p output
	../bin/halRandGen --preset small --seed 0 --testRand --format hdf5 output/small.hdf5.hal

../bin/halRandGen:
	cd ../randgen && ${MAKE}

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

