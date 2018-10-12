binDir =${rootDir}/bin
libDir = ${rootDir}/lib
objDir = ${rootDir}/objs

libHal = ${libDir}/libHal.a
libHalStats = ${libDir}/libHalStats.a
libHalChain = ${libDir}/libHalChain.a
libHalMutations = ${libDir}/libHalMutations.a
libHalLiftover = ${libDir}/libHalLiftover.a
libHalLod = ${libDir}/libHalLod.a
libHalMaf = ${libDir}/libHalMaf.a

inclSpec = -I${rootDir}/api/inc -Iimpl -Iinc

#Modify this variable to set the location of sonLib
#(sonlib is used only for cuTest at this potin)
sonLibRootDir ?= ${rootDir}/../sonLib
sonLibDir = ${sonLibRootDir}/lib

# update PYTHONPATH for tests
export PYTHONPATH := $(abspath ${rootDir}/..):${PYTHONPATH}

export PATH := $(abspath ${rootDir}/bin):${PATH}

.SECONDARY: 

include  ${sonLibRootDir}/include.mk

dataSetsDir=/Users/hickey/Documents/Devel/genomes/datasets

#
# The -D_GLIBCXX_USE_CXX11_ABI=1 flag prevents errors with
#   std::__cxx11::basic_string vs std::basic_string
# when linking with HDF5 and modern C++ compilers (well, GCC 7.3)
#
cflags += -I${sonLibDir} -fPIC
cppflags += -I${sonLibDir} -fPIC -D_GLIBCXX_USE_CXX11_ABI=1 -std=c++11 -Wno-sign-compare

basicLibs = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a
basicLibsDependencies = ${basicLibs}

# hdf5 compilation is done through its wrappers.
# we can speficy our own (sonlib) compilers with these variables:
cpp = h5c++ ${h5prefix}
cxx = h5cc ${h5prefix}

# add compiler flag and kent paths if udc is enabled
# relies on KENTSRC containing path to top level kent/ dir
# and MACHTYPE being specified
ifdef ENABLE_UDC
#  Find samtabix as in kent/src/inc/common.mk:
	ifeq (${SAMTABIXDIR},)
		SAMTABIXDIR = /hive/data/outside/samtabix/${MACHTYPE}
	endif

	cppflags += -DENABLE_UDC -I${KENTSRC}/src/inc -pthread
	cflags += -I${KENTSRC}/src/inc -pthread
	basicLibs += ${KENTSRC}/src/lib/${MACHTYPE}/jkweb.a  ${SAMTABIXDIR}/libsamtabix.a -lssl -lcrypto
endif


#
# phyloP support
#
phyloPcppflags = 
phyloPlibs = 

ifdef ENABLE_PHYLOP

ifndef TARGETOS
  TARGETOS := $(shell uname -s)
endif

#	Local Linux install (phast and clapack sister dirs to hal/)
#	(note CLAPACKPATH not needed in Mac)
	PHAST=../../phast
	CLAPACKPATH=../../clapack

ifeq ($(TARGETOS), Darwin)
	phyloPcppflags += -DENABLE_PHYLOP -I${PHAST}/include -I${PHAST}/src/lib/pcre -DVECLIB
	phyloPlibs += -L${PHAST}/lib -lphast -lc -framework Accelerate
else
	F2CPATH=${CLAPACKPATH}/F2CLIBS
	phyloPcppflags += -DENABLE_PHYLOP -I${PHAST}/include -I${PHAST}/src/lib/pcre -I${CLAPACKPATH}/INCLUDE -I${F2CPATH}
	phyloPlibs += -L${PHAST}/lib -lphast -L${CLAPACKPATH} -L${F2CPATH} -llapack -ltmg -lblaswr -lf2c 
endif

endif

# test includes and libs uses buy several modules
halApiTestIncl = ${rootDir}/api/tests
halApiTestSupportLibs = ${objDir}/api/tests/halAlignmentTest.o ${objDir}/api/tests/halBottomSegmentTest.o \
    ${objDir}/api/tests/halTopSegmentTest.o ${objDir}/api/tests/halAlignmentInstanceTest.o \
    ${objDir}/api/tests/halRandomData.o
