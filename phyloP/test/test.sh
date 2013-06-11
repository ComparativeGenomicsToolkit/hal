hal2maf blanchette.hal blanchette.maf --refGenome HUMAN
halPhyloP blanchette.mod blanchette.hal HUMAN  > fromHal
phyloP -i MAF --method LRT --mode CONACC --wig-scores blanchette.mod blanchette.maf > fromMaf
