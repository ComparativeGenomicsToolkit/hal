rootDir = .
include include.mk

modules = api stats randgen validate mutations fasta alignmentDepth liftover lod maf blockViz extract analysis phyloP modify assemblyHub synteny paf


.PHONY: all libs %.libs progs %.progs clean %.clean doxy %.doxy

all : libs progs

libs: ${modules:%=%.libs}
%.libs:
	cd $* && ${MAKE} libs

progs: ${modules:%=%.progs}
%.progs: libs
	cd $* && ${MAKE} progs

clean: ${modules:%=%.clean}
	rm -f hal
	rm -rf lib bin objs
	rm -f *.pyc */*.pyc */*/*.pyc
	rm -rf __pycache__ */__pycache__ */*/__pycache__

%.clean:
	cd $* && ${MAKE} clean

# create symbolic links for test so that python packages work without assuming name of
# directory, but then remove, as it can cause grief for naive copy programs
test:
	rm -f hal
	ln -sf . hal
	${MAKE} doTests
	rm -f hal


doTests: ${modules:%=%.test}

%.test: all
	cd $* && ${MAKE} test


doxy : ${modules:%=doxy.%}

doxy.%:
	cd api && ${MAKE} doxy

etags:
	etags $$(find . -name '*.h' -o -name '*.cpp')
