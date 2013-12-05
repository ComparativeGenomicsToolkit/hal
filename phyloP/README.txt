After unpacking this file, move the directory to the hal source code directory:
tar -zxf phyloP.tar.gz
mv phyloP $HALPATH/

Need to install phast, which requires clapack:
 wget http://www.netlib.org/clapack/clapack.tgz
 tar -xvzf clapack.tgz
 cd CLAPACK-3.2.1
 cp make.inc.example make.inc && make f2clib && make blaslib && make lib
 cd ..
 svn co http://compgen.bscb.cornell.edu/svnrepo/phast/trunk phast/
 cd phast/src
 make CLAPACKPATH=/path/to/CLAPACK-3.2.1

in $HALPATH/phyloP:
edit Makefile to point to PHAST and CLAPACKPATH. (F2CPATH should not need to be changed.)
type 'make', should create binary halPhyloP

