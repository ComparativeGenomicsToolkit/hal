#!/bin/bash                                                                    

binSize=100000
reps=1000
ref=reference
tgts=EscherichiaColiWUid162011,EscherichiaColiKo11flUid162099,EscherichiaColiKo11flUid52593,ShigellaSonnei53gUid84383
refChr=referencerefChr1
len=9808498
lodPath=http://hgwdev.cse.ucsc.edu/~nknguyen/ecoli/hub/pangenome-dups/lod.txt
halPath=http://hgwdev.cse.ucsc.edu/~nknguyen/ecoli/hub/pangenome-dups/out.hal
udcPath=/hive/users/hickey/udcTemp
outPath=/hive/users/hickey/timing

./blockVizBenchmark.py ${lodPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} > ${outPath}/ecoli.csv

./blockVizBenchmark.py ${lodPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} --zapUdc > ${outPath}/ecoli_noudc.csv

#./blockVizBenchmark.py ${halPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} > ${outPath}/ecoli_nolod.csv