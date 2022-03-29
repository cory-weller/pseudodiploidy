singularity exec -H $PWD:/home src/R.sif Rscript \
    src/buildCountMatrix.R \
    data/input/combined_featurecounts.csv \
    data/input/samples.csv \
    data/processed/featurecounts.mat.RDS