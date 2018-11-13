#!/bin/bash                                                                    

binSize=1000000
reps=1000
ref=reference
tgts=EscherichiaColiWUid162011,EscherichiaColiKo11flUid162099,EscherichiaColiKo11flUid52593,ShigellaSonnei53gUid84383
refChr=referencerefChr1
len=9808498
lodPath=http://hgwdev.cse.ucsc.edu/~nknguyen/ecoli/hub/pangenome-dups/lod.txt
halPath=http://hgwdev.cse.ucsc.edu/~nknguyen/ecoli/hub/pangenome-dups/out.hal
udcPath=/hive/users/hickey/udcTemp
outPath=/hive/users/hickey/timing
seed=2323
python=/cluster/home/jcarmstr/Python-2.7/python

${python} ./blockVizBenchmark.py ${lodPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} --seed ${seed} > ${outPath}/ecoli.csv

${python} ./blockVizBenchmark.py ${lodPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} --zapUdc --seed ${seed} > ${outPath}/ecoli_noudc.csv

${python} ./blockVizBenchmark.py ${halPath} ${ref} ${refChr} ${len} ${tgts} --udc ${udcPath} --reps ${reps} --binSize ${binSize} --seed ${seed} > ${outPath}/ecoli_nolod.csv