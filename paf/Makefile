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

hal2pafSmallMMapTest: tests/output/small.mmap1.0.hal
	../bin/hal2paf tests/output/small.mmap1.0.hal --onlySequenceNames | gzip > tests/output/$@.paf.gz
	zdiff tests/expected/$@.paf.gz tests/output/$@.paf.gz

hal2pafMouseRatTest:
	../bin/hal2paf tests/input/mr.hal | gzip > tests/output/$@.paf.gz
	zdiff tests/expected/$@.paf.gz tests/output/$@.paf.gz

tests/output/small.mmap1.0.hal:
	mkdir -p tests/output
	bunzip2 -dc ../extract/tests/input/small.mmap1.0.hal.bz2 > tests/output/small.mmap1.0.hal

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}



# Local Variables:
# mode: makefile-gmake
# End:

