#!/usr/bin/env bash

chromosome=${1}
start=${2}
stop=${3}

vcf="data/external/chromosome${chromosome}.vcf.gz"

singularity exec --bind $PWD src/pseudodiploidy.sif tabix ${vcf} chromosome${chromosome}:${start}-${stop} 2>/dev/null | grep -F "chromosome${chromosome}"
