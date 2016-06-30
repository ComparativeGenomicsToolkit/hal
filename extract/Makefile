rootPath = ../
include ../include.mk

libTestSources = $(wildcard tests/*.cpp)
libTestHeaders = $(wildcard tests/*.h)
libTestsCommon = ${rootPath}/api/tests/halAlignmentTest.cpp ${rootPath}/api/tests/halAlignmentInstanceTest.cpp
libTestsCommonHeaders = ${rootPath}/api/tests/halAlignmentTest.h ${rootPath}/api/tests/halAlignmentInstanceTest.h ${rootPath}/api/tests/allTests.h

all : ${binPath}/halExtract ${binPath}/halAlignedExtract ${binPath}/halMaskExtract ${binPath}/hal4dExtract ${binPath}/hal4dExtractTest ${binPath}/halSingleCopyRegionsExtract

clean : 
	rm -f ${binPath}/halExtract ${binPath}/halAlignedExtract ${binPath}/halMaskExtract ${binPath}/hal4dExtract

${binPath}/halExtract : impl/halExtract.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halExtract impl/halExtract.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halAlignedExtract :impl/halAlignedExtract.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halAlignedExtract impl/halAlignedExtract.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/halMaskExtract : impl/halMaskExtractMain.cpp  impl/halMaskExtractor.cpp inc/halMaskExtractor.h ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halMaskExtract impl/halMaskExtractMain.cpp  impl/halMaskExtractor.cpp ${libPath}/halLib.a ${basicLibs}

${binPath}/hal4dExtract : impl/hal4dExtractMain.cpp  impl/hal4dExtract.cpp inc/hal4dExtract.h ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -o ${binPath}/hal4dExtract impl/hal4dExtractMain.cpp  impl/hal4dExtract.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

${binPath}/halSingleCopyRegionsExtract : impl/halSingleCopyRegionsExtract.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -o ${binPath}/halSingleCopyRegionsExtract impl/halSingleCopyRegionsExtract.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}

${binPath}/hal4dExtractTest : impl/hal4dExtract.cpp tests/hal4dExtractTest.cpp ${libTestSources} ${libTestHeaders} ${libTestsCommon} ${libTestsHeadersCommon} ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I tests -I ../api/tests -o ${binPath}/hal4dExtractTest impl/hal4dExtract.cpp ${libTestsCommon} tests/hal4dExtractTest.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}
