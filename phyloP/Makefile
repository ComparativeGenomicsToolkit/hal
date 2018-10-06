rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/phyloP

halPhyloP_srcs = impl/halPhyloPMain.cpp impl/halPhyloPBed.cpp impl/halPhyloP.cpp
halPhyloP_objs = ${halPhyloP_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${halPhyloP_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs =  ${binDir}/halPhyloP ${binDir}/halPhyloPTrain.py ${binDir}/halPhyloPMP.py ${binDir}/halTreePhyloP.py

ifdef ENABLE_PHYLOP
all: progs
libs:
progs: ${progs}
else
all:
libs:
progs:
endif

clean: 
	rm -f ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

