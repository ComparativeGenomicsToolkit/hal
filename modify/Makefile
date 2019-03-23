rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/modify

markAncestors_srcs = markAncestors.cpp
markAncestors_objs = ${markAncestors_srcs:%.cpp=${modObjDir}/%.o}
renameFile_srcs = renameFile.cpp
renameFile_objs = ${renameFile_srcs:%.cpp=${modObjDir}/%.o}
halRemoveGenome_srcs = halRemoveGenome.cpp
halRemoveGenome_objs = ${halRemoveGenome_srcs:%.cpp=${modObjDir}/%.o} ${markAncestors_objs}
halAddToBranch_srcs = halAddToBranch.cpp
halAddToBranch_objs = ${halAddToBranch_srcs:%.cpp=${modObjDir}/%.o} ${markAncestors_objs}
halReplaceGenome_srcs = halReplaceGenome.cpp
halReplaceGenome_objs = ${halReplaceGenome_srcs:%.cpp=${modObjDir}/%.o} ${markAncestors_objs}
halReplaceSequence_srcs = halReplaceSequence.cpp
halReplaceSequence_objs = ${halReplaceSequence_srcs:%.cpp=${modObjDir}/%.o}
halAppendSubtree_srcs = halAppendSubtree.cpp
halAppendSubtree_objs = ${halAppendSubtree_srcs:%.cpp=${modObjDir}/%.o} ${markAncestors_objs}
findRegionsExclusivelyInGroup_srcs = findRegionsExclusivelyInGroup.cpp
findRegionsExclusivelyInGroup_objs = ${findRegionsExclusivelyInGroup_srcs:%.cpp=${modObjDir}/%.o}
halUpdateBranchLengths_srcs = halUpdateBranchLengths.cpp
halUpdateBranchLengths_objs = ${halUpdateBranchLengths_srcs:%.cpp=${modObjDir}/%.o}
halWriteNucleotides_srcs = halWriteNucleotides.cpp
halWriteNucleotides_objs = ${halWriteNucleotides_srcs:%.cpp=${modObjDir}/%.o}
halSetMetadata_srcs = halSetMetadata.cpp
halSetMetadata_objs = ${halSetMetadata_srcs:%.cpp=${modObjDir}/%.o}
halRenameGenomes_srcs = halRenameGenomes.cpp
halRenameGenomes_objs = ${halRenameGenomes_srcs:%.cpp=${modObjDir}/%.o} ${renameFile_objs}
halRenameSequences_srcs = halRenameSequences.cpp
halRenameSequences_objs = ${halRenameSequences_srcs:%.cpp=${modObjDir}/%.o} ${renameFile_objs}
ancestorsML_srcs = ancestorsML.cpp
ancestorsML_objs = ${ancestorsML_srcs:%.cpp=${modObjDir}/%.o}
halAncestorsMLTest_srcs = halAncestorsMLTest.cpp
halAncestorsMLTest_objs = ${halAncestorsMLTest_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${markAncestors_srcs} ${renameFile_srcs} ${halRemoveGenome_srcs} ${halAddToBranch_srcs} \
    ${halReplaceGenome_srcs} ${halAppendSubtree_srcs} \
    ${findRegionsExclusivelyInGroup_srcs} ${halUpdateBranchLengths_srcs} \
    ${halWriteNucleotides_srcs} ${halSetMetadata_srcs} ${halRenameGenomes_srcs} \
    ${halRenameSequences_srcs} ${ancestorsML_srcs} ${halAncestorsMLTest_srcs} ${halReplaceSequence_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halRemoveGenome ${binDir}/halAddToBranch ${binDir}/halReplaceGenome ${binDir}/halAppendSubtree ${binDir}/findRegionsExclusivelyInGroup ${binDir}/halUpdateBranchLengths ${binDir}/halWriteNucleotides ${binDir}/halSetMetadata ${binDir}/halRenameGenomes ${binDir}/halRenameSequences ${binDir}/halReplaceSequence

inclSpec += -I${rootDir}/liftover/inc
otherLibs += ${libHalLiftover}
ifdef ENABLE_PHYLOP
inclSpec += ${phyloPcppflags}
otherLibs += ${phyloPlibs}
phast_progs = ${binDir}/ancestorsML ${binDir}/halAncestorsMLTest ${binDir}/ancestorsMLMP.py
endif

ifdef ENABLE_PHYLOP
all: progs
libs:
progs: ${progs} ${phast_progs}
else
all: progs
libs:
progs: ${progs}
endif

clean : 
	rm -f ${objs} ${progs} ${phast_progs} ${depends}
test:

include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

