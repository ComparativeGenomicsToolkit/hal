rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/api

libHal_srcs = $(wildcard impl/*.cpp) $(wildcard hdf5_impl/*.cpp) $(wildcard mmap_impl/*.cpp)
ifdef ENABLE_UDC
  libHal_c_srcs = impl/udc2.c
endif
libHal_objs = ${libHal_srcs:%.cpp=${modObjDir}/%.o} ${libHal_c_srcs:%.c=${modObjDir}/%.o}
halApiTests_srcs = $(wildcard tests/*.cpp)
halApiTests_objs = ${halApiTests_srcs:%.cpp=${modObjDir}/%.o}
halHdf5Tests_srcs = $(wildcard hdf5_tests/*.cpp)
halHdf5Tests_objs = ${halHdf5Tests_srcs:%.cpp=${modObjDir}/%.o}
ifdef ENABLE_UDC
   udc2Tests_srcs = $(wildcard tests/udc2Test.c)
   udc2Tests_objs = ${udc2Tests_srcs:%.c=${modObjDir}/%.o}
endif
srcs = ${libHal_srcs} ${halApiTests_srcs} ${halHdf5Tests_srcs}
c_srcs = ${udc2Tests_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o} ${c_srcs:%.c=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend} ${c_srcs:%.c=%.depend}

progs =  ${binDir}/halApiTests ${binDir}/halHdf5Tests
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

halApiTests: hdf5.halApiTests mmap.halApiTests

%.halApiTests:
	${binDir}/halApiTests $*

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
