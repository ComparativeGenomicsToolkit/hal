rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/api

libHal_srcs = $(wildcard impl/*.cpp) $(wildcard hdf5_impl/*.cpp) $(wildcard mmap_impl/*.cpp)
ifdef ENABLE_UDC
  libHal_c_srcs = impl/udc2.c
endif
libHal_objs = ${libHal_srcs:%.cpp=${modObjDir}/%.o} ${libHal_c_srcs:%.c=${modObjDir}/%.o}

halHdf5Tests_srcs = $(wildcard hdf5_tests/*.cpp)
halHdf5Tests_objs = ${halHdf5Tests_srcs:%.cpp=${modObjDir}/%.o}
halHdf5Tests_progs = ${binDir}/halHdf5Tests

halApiTests_srcs = $(wildcard tests/*.cpp)
halApiTests_objs = ${halApiTests_srcs:%.cpp=${modObjDir}/%.o}

halApiTest_names = halAlignmentTreesTest \
	halBottomSegmentTest \
	halColumnIteratorTest \
	halGappedSegmentIteratorTest \
	halGenomeTest \
	halMappedSegmentTest \
	halMetaDataTest \
	halRearrangementTest \
	halSequenceTest \
	halTopSegmentTest \
	halValidateTest
halApiTest_progs = ${halApiTest_names:%=${binDir}/%}

# make magic to generate the variables containing the objects for the link rule.
# for each prog name this generates a _objs variable (e.g. halValidateTest_objs)
$(foreach prog,${halApiTest_names},$(eval ${prog}_objs = ${modObjDir}/tests/${prog}.o ${halApiTestSupportLibs}))

ifdef ENABLE_UDC
   udc2Tests_srcs = $(wildcard tests/udc2Test.c)
   udc2Tests_objs = ${udc2Tests_srcs:%.c=${modObjDir}/%.o}
endif
srcs = ${libHal_srcs} ${halApiTests_srcs} ${halHdf5Tests_srcs}
c_srcs = ${udc2Tests_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o} ${c_srcs:%.c=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend} ${c_srcs:%.c=%.depend}

progs = ${halHdf5Tests_progs} ${halApiTest_progs}
inclSpec += -Ihdf5_impl -Immap_impl
ifdef ENABLE_UDC
   # FIXME: standarize var names
   inclSpec += -I${KENTSRC}/inc -I${KENTSRC}/htslib
   progs +=  ${binDir}/udc2Tests
endif


all : libs progs
libs: ${libHal} ${halApiTestSupportLibs} 
progs: ${progs}

clean : 
	rm -f ${libHal} ${objs} ${progs} ${depends}

test: halHdf5Tests halApiTests

halHdf5Tests: all
	${binDir}/halHdf5Tests


halApiTests: hdf5.halApiTestsStorage mmap.halApiTestsStorage

%.halApiTestsStorage:
	${MAKE} runHalApiTest halStorageFormat=$*

runHalApiTest: ${halApiTest_names:%=%.runHalApiTest}

%.runHalApiTest:
	${binDir}/$* ${halStorageFormat}

doxy :
	doxygen doc/doxy.cfg

etags:
	(cd .. && ${MAKE} etags}

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:
