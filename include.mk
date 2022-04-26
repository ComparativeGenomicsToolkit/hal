# if include.local.mk exists, include it first to set various options
# it shuld not be checked in
includeLocal = ${rootDir}/include.local.mk
ifneq ($(wildcard ${includeLocal}),)
   include ${includeLocal}
endif

binDir =${rootDir}/bin
libDir = ${rootDir}/lib
objDir = ${rootDir}/objs

libHal = ${libDir}/libHal.a
libHalStats = ${libDir}/libHalStats.a
libHalBlockViz = ${libDir}/libHalBlockViz.a
libHalMutations = ${libDir}/libHalMutations.a
libHalLiftover = ${libDir}/libHalLiftover.a
libHalLod = ${libDir}/libHalLod.a
libHalMaf = ${libDir}/libHalMaf.a

inclSpec += -I${rootDir}/api/inc -Iimpl -Iinc -I${rootDir}/liftover/inc

#Modify this variable to set the location of sonLib
#(sonlib is used only for cuTest at this potin)
sonLibRootDir ?= ${rootDir}/../sonLib
sonLibDir = ${sonLibRootDir}/lib

# update PYTHONPATH and PATH for tests, symlink 'hal' allows checkout directory to have
# any name
export PYTHONPATH := $(abspath ${sonLibRootDir}/src):$(abspath ${rootDir}):${PYTHONPATH}

export PATH := $(abspath ${rootDir}/bin):${PATH}

.SECONDARY: 

include  ${sonLibRootDir}/include.mk

##
# For GCC, the C++ 11 aplication binary interface must match the version
# that HDF5 uses.  Change the ABI version can address errors about undefined functions that vary
# by string type, such as 
#   std::__cxx11::basic_string vs std::basic_string
#
# Specify one of:
#   CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=0
# or
#   CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=1
#
# in include.local.mk or the environment will change the API.
# You must do a make clean if you change this variable.

ifeq (${CXX_ABI_DEF},)
    CXX_ABI_DEF = -D_GLIBCXX_USE_CXX11_ABI=1
endif

CFLAGS += -I${sonLibDir}
CXXFLAGS += -I${sonLibDir} ${CXX_ABI_DEF} -std=c++11 -Wno-sign-compare

LDLIBS += ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a
LIBDEPENDS += ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a

# hdf5 compilation is done through its wrappers.  See README.md for discussion of
# h5prefix
CXX = h5c++ ${h5prefix}
CC = h5cc ${h5prefix}

#
# phyloP support
#
phyloPCXXFLAGS = 
phyloPlibs = 

ifdef ENABLE_PHYLOP

ifndef TARGETOS
  TARGETOS := $(shell uname -s)
endif

#  Defaults to local Linux install (phast and clapack sister dirs to hal/)
# (note CLAPACKPATH not needed in Mac)
ifeq (${PHAST},)
    PHAST=${rootDir}/../phast
endif
ifeq (${CLAPACKPATH},)
    CLAPACKPATH=${rootDir}/../clapack
endif

# pointing at both ${PHAST}/include/phast and ${PHAST}/include allows compiling
# with v1.6 or v1.5.  However, only do the includes where needed, to avoid
# include file name conflicts.
ifeq ($(TARGETOS), Darwin)
    PHASTCXXFLAGS += -DENABLE_PHYLOP -I${PHAST}/include/phast -I${PHAST}/include -I${PHAST}/src/lib/pcre -DVECLIB
    LDLIBS += -L${PHAST}/lib -lphast -lc -framework Accelerate
else
    F2CPATH=${CLAPACKPATH}/F2CLIBS
    PHASTCXXFLAGS += -DENABLE_PHYLOP -I${PHAST}/include/phast -I${PHAST}/include -I${PHAST}/src/lib/pcre -I${CLAPACKPATH}/INCLUDE -I${F2CPATH}
    LDLIBS += -L${PHAST}/lib -lphast -L${CLAPACKPATH} -L${F2CPATH} -llapack -ltmg -lblaswr -lf2c 
endif

endif

# add compiler flag and kent paths if udc is enabled
# relies on KENTSRC containing path to top level kent/src dir
# and MACHTYPE being specified.
# This MUST follow PHAST defs, as they both have a gff.h
ifdef ENABLE_UDC
    #  Find htslib as in kent/src/inc/common.mk:
    MACHTYPE = x86_64
    CXXFLAGS += -DENABLE_UDC
    CFLAGS += -DENABLE_UDC
    UDCCXXFLAGS += -I${KENTSRC}/inc -I${KENTSRC}/htslib -pthread
    UDCCFLAGS += -Wall -Werror -std=c99 -I${KENTSRC}/inc -I${KENTSRC}/htslib
    LDLIBS += ${KENTSRC}/lib/${MACHTYPE}/jkweb.a  ${KENTSRC}/htslib/libhts.a -lcurl -lssl -lcrypto -pthread
endif



# test includes and libs uses buy several modules
halApiTestIncl = ${rootDir}/api/tests
halApiTestSupportLibs = ${objDir}/api/tests/halApiTestSupport.o ${objDir}/api/tests/halRandomData.o

LDLIBS += ${LIBS}
