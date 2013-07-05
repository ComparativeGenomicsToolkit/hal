hal2maf blanchette.hal blanchette.maf --refGenome HUMAN
halPhyloP blanchette.hal HUMAN blanchette.mod fromHal
phyloP -i MAF --method LRT --mode CONACC --wig-scores blanchette.mod blanchette.maf > fromMaf
