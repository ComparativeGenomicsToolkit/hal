rootDir = ..
include ${rootDir}/include.mk

progs = ${binDir}/hal2assemblyHub.py

all: libs progs
libs:
progs: ${progs}
clean:
	rm -f ${progs}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

