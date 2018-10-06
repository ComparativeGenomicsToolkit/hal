rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/api

progs = ${binDir}/halTreeNIConservation.py  ${binDir}/halTreeNIBackground.py ${binDir}/halTreeNITurnover.py

all: progs
libs:
progs:  ${progs}
clean: 
	rm -f ${progs}
test:

${binDir}/%.py : neutralIndel/%.py
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod u+x,a-w $@

# Local Variables:
# mode: makefile-gmake
# End:

