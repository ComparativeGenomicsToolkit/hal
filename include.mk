#Location of sonLib
binPath=${rootPath}bin/
libPath=${rootPath}lib/

#Modify this variable to set the location of sonLib
#(sonlib is used only for cuTest at this potin)
sonLibRootPath=${rootPath}/../sonLib
sonLibPath=${sonLibRootPath}/lib

include  ${sonLibRootPath}/include.mk

dataSetsPath=/Users/hickey/Documents/Devel/genomes/datasets

cflags += -I${sonLibPath}
cppflags += -I${sonLibPath}

basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a
basicLibsDependencies = ${sonLibPath}/cuTest.a 

# hdf5 compilation is done through its wrappers.
# we can speficy our own (sonlib) compilers with these variables:
HDF5_CXX = ${cpp}
HDF5_CXXLINKER = ${cpp}
HDF5_CC = ${cxx}
HDF5_CCLINKER = ${cxx} 
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
