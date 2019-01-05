#!/bin/bash
# Example invocation of the Regulatory Locus Intersection (RELI) tool using
# sample data in ./example (European ancestry)

RELIBIN=${RELI_BINARY:-../RELI}
DATADIR=${RELI_DATA_DIR:-../data}
OUTPUTDIR=${RELI_OUTPUT_DIR:-../output}

OPTIONS=(
    -snp SLE_EU.snp
    -ld SLE_EU.ld
    -index "$DATADIR/ChIPseq.index"
    -data "$DATADIR/ChIP-seq"
    -target hg19_0302
    -build "$DATADIR/GenomeBuild/hg19.txt"
    -null "$DATADIR/Null/CommonSNP_MAFmatch"
    -dbsnp "$DATADIR/SNPtable/SNPtable"
    -out "$OUTPUTDIR"
    -match
    -rep 2000
    -corr 1544
    -phenotype Systemic_Lupus_Erythematosus
    -ancestry EU
)

# invoke RELI; double quotes are essential here in case any of the options
# above contain embedded spaces
"$RELIBIN" "${OPTIONS[@]}"

