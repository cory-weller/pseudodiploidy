#!/usr/bin/env bash

strainsFile=${1}
chromosome=${2}
start=${3}
stop=${4}
replicates=${5}

singularity exec --bind $PWD src/pseudodiploidy.sif Rscript src/query-region.R ${strainsFile} ${chromosome} ${start} ${stop} ${replicates}
