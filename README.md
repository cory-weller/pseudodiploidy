# Yeast Pseudodiploidy

## Singularity
### Install Singularity
```
bash src/install-singularity.sh
```

### Build Singularity image
```
singularity build --fakeroot src/R.sif src/R.def
```

### Upload Singularity Image 
`R.sif` uploaded to OneDrive, which can be downloaded with:
```
wget -O R.sif "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119114&authkey=AH63WVNezJbyWqM" 
```

## Retrieve data
```
bash src/get-data.sh
```

## Transform data
```
# Build count matrix for DESeq
singularity exec --bind ${PWD} src/R.sif Rscript \
    src/build-count-matrix.R \
    data/input/combined-featurecounts.csv \
    data/input/samples.csv \
    data/processed/featurecounts-matrix.RDS

# Build TPM table for other analyses
singularity exec --bind ${PWD} src/R.sif Rscript \
    src/build-TPM-table.R \
    data/input/combined-featurecounts.csv \
    data/input/samples.csv \
    data/processed/TPM.txt.gz

# Build DESeq Data Set (DDS) Object
singularity exec --bind ${PWD} src/R.sif Rscript \
    src/build-DDS.R \
    data/processed/featurecounts-matrix.RDS \
    data/input/samples.csv \
    data/processed/DDS.RDS

# Run differential expression contrasts
# using src/run-contrast.R
# (run manually in R for now)
# by starting R session:
./R.sh

```

# Run exploratory analyses
```
singularity exec --bind ${PWD} src/R.sif R
```
