rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/

hal2fasta_srcs = hal2fasta.cpp 
hal2fasta_objs = ${hal2fasta_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${hal2fasta_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/hal2fasta

all: progs
libs:
progs: ${progs}

clean: 
	rm -f  ${objs} ${progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

