rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/api

libHal_srcs = $(wildcard impl/*.cpp) $(wildcard hdf5_impl/*.cpp) $(wildcard mmap_impl/*.cpp)
libHal_objs =  ${libHal_srcs:%.cpp=${modObjDir}/%.o}
halApiTests_srcs = $(wildcard tests/*.cpp)
halApiTests_objs = ${halApiTests_srcs:%.cpp=${modObjDir}/%.o}
halHdf5Tests_srcs = $(wildcard hdf5_tests/*.cpp)
halHdf5Tests_objs = ${halHdf5Tests_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHal_srcs} ${halApiTests_srcs} ${halHdf5Tests_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
inclSpec += -Ihdf5_impl -Immap_impl
depends = ${srcs:%.cpp=%.depend}

progs =  ${binDir}/halApiTests ${binDir}/halHdf5Tests

all : libs progs
libs: ${libHal} ${halApiTestSupportLibs} 
progs: ${progs}

clean : 
	rm -f ${libHal} ${objs} ${progs} ${depends}

test: halHdf5Tests halApiTests

halHdf5Tests:
	${binDir}/halHdf5Tests
halApiTests:
	${binDir}/halApiTests

doxy :
	doxygen doc/doxy.cfg

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:
