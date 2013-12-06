rootPath = ../
include ../include.mk

all: ${binPath}/hal2assemblyHub.py

clean:
	rm -f ${binPath}/hal2assemblyHub.py

#${binPath}/hal2assemblyHub.py : hal2assemblyHub.py bedCommon.py gcPercentTrack.py prepareHubFiles.py snakeTrack.py alignabilityTrack.py bedTrack.py groupExclusiveRegions.py prepareLodFiles.py treeCommon.py assemblyHubCommon.py conservationTrack.py rmskTrack.py wigTrack.py docs/conservationDocs.py docs/hubCentralDocs.py docs/repeatMaskerDocs.py docs/alignabilityDocs.py docs/gcPercentDocs.py docs/makeDocs.py
${binPath}/hal2assemblyHub.py : hal2assemblyHub.py
	cp hal2assemblyHub.py ${binPath}/hal2assemblyHub.py
	chmod +x ${binPath}/hal2assemblyHub.py


