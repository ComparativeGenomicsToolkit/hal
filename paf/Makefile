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

test: hal2pafSmallMMapTest hal2pafMouseRatTest

hal2pafSmallMMapTest: tests/output/small.mmap1.0.hal tests/output/hal2pafSmallMMapTest.paf.baseline
	../bin/hal2paf tests/output/small.mmap1.0.hal --onlySequenceNames > tests/output/$@.paf
	diff tests/output/$@.paf tests/output/hal2pafSmallMMapTest.paf.baseline

hal2pafMouseRatTest: tests/output/hal2pafMouseRatTest.paf.baseline
	../bin/hal2paf tests/input/mr.hal > tests/output/$@.paf
	diff tests/output/$@.paf tests/output/hal2pafMouseRatTest.paf.baseline

tests/output/small.mmap1.0.hal: output
	bunzip2 -dc ../extract/tests/input/small.mmap1.0.hal.bz2 > tests/output/small.mmap1.0.hal

tests/output/hal2pafSmallMMapTest.paf.baseline: output
	gzip -dc tests/expected/hal2pafSmallMMapTest.paf.gz > tests/output/hal2pafSmallMMapTest.paf.baseline

tests/output/hal2pafMouseRatTest.paf.baseline: output
	gzip -dc tests/expected/hal2pafMouseRatTest.paf.gz > tests/output/hal2pafMouseRatTest.paf.baseline

output:
	mkdir -p tests/output

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}



# Local Variables:
# mode: makefile-gmake
# End:

