rootPath = ../
include ../include.mk

all: ${binPath}/hal2assemblyHub.py

clean:
	rm -f ${binPath}/hal2assemblyHub.py

${binPath}/hal2assemblyHub.py : hal2assemblyHub.py
	cp hal2assemblyHub.py ${binPath}/hal2assemblyHub.py
	chmod +x ${binPath}/hal2assemblyHub.py
