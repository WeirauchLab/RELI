#!/bin/bash
# Example invocation of the Regulatory Locus Intersection (RELI) tool

../RELI \
  -snp SLE_EU.snp  \
  -ld SLE_EU.ld \
  -index ../data/ChIPseq.index \
  -data ../data/ChIP-seq \
  -target hg19_0302 \
  -build ../data/GenomeBuild/hg19.txt \
  -null ../data/Null/CommonSNP_MAFmatch \
  -dbsnp ../data/SNPtable/SNPtable \
  -out Output   \
  -match \
  -rep 2000 \
  -corr 1544 \
  -phenotype Systemic_Lupus_Erythematosus \
  -ancestry EU
