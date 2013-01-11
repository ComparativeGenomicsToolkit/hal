rootPath = ../
include ../include.mk

all : ${binPath}/halTreeNIConservation.py  ${binPath}/halTreeNIBackground.py ${binPath}/halTreeNITurnover.py

clean : 
	rm -f  ${binPath}/halTreeNIConservation.py  ${binPath}/halTreeNIBackground.py ${binPath}/halTreeNITurnover.py

${binPath}/halTreeNIConservation.py : neutralIndel/halTreeNIConservation.py
	cp neutralIndel/halTreeNIConservation.py ${binPath}
	chmod u+x ${binPath}/halTreeNIConservation.py

${binPath}/halTreeNIBackground.py : neutralIndel/halTreeNIBackground.py
	cp neutralIndel/halTreeNIBackground.py ${binPath}
	chmod u+x ${binPath}/halTreeNIBackground.py

${binPath}/halTreeNITurnover.py : neutralIndel/halTreeNITurnover.py
	cp neutralIndel/halTreeNITurnover.py ${binPath}
	chmod u+x ${binPath}/halTreeNITurnover.py
