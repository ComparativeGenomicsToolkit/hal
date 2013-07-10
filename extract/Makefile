rootPath = ../
include ../include.mk

all : ${binPath}/halExtract ${binPath}/halAlignedExtract ${binPath}/halMaskExtract ${binPath}/hal4dExtract

clean : 
	rm -f ${binPath}/halExtract ${binPath}/halAlignedExtract ${binPath}/halMaskExtract ${binPath}/hal4dExtract

${binPath}/halExtract : impl/halExtract.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halExtract impl/halExtract.cpp ${libPath}/halLib.a ${basicLibs}
	rm *.o

${binPath}/halAlignedExtract :impl/halAlignedExtract.cpp ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I ${rootPath}/api/tests -o ${binPath}/halAlignedExtract impl/halAlignedExtract.cpp ${libPath}/halLib.a ${basicLibs}
	rm *.o

${binPath}/halMaskExtract : impl/halMaskExtractMain.cpp  impl/halMaskExtractor.cpp inc/halMaskExtractor.h ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/halMaskExtract impl/halMaskExtractMain.cpp  impl/halMaskExtractor.cpp ${libPath}/halLib.a ${basicLibs}
	rm *.o

${binPath}/hal4dExtract : impl/hal4dExtractMain.cpp  impl/hal4dExtract.cpp inc/hal4dExtract.h ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibsDependencies}
	${cpp} ${cppflags} -I inc -I impl -I ${libPath} -I impl -I tests -o ${binPath}/hal4dExtract impl/hal4dExtractMain.cpp  impl/hal4dExtract.cpp ${libPath}/halLiftover.a ${libPath}/halLib.a ${basicLibs}
	rm *.o
