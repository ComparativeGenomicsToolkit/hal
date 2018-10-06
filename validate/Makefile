rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/validate

halValidate_srcs = halValidateMain.cpp
halValidate_objs = ${halValidate_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${halValidate_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halValidate

inclSpec += -I${rootDir}/stats/inc/


all: progs
libs:
progs: ${progs}

clean: 
	rm -f ${objs} ${progs} ${depends}
test:


include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

