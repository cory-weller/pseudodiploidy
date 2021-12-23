# Yeast Pseudodiploidy

## Singularity
### Install Singularity
```
bash src/install_singularity.sh
```

### Build Singularity image
```
singularity build --fakeroot src/R.sif src/r.def
```

### Upload dd Singularity Image 
`R.sif` uploaded to OneDrive, which can be downloaded with:
```
wget -O R.sif "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2119110&authkey=AH7s36_oIlV8uIU"
```

## Retrieve data
```
bash src/getData.sh
```

## Transform data
```
# Build count matrix for DESeq
singularity exec -H $PWD:/home src/R.sif Rscript \
    src/buildCountMatrix.R \
    data/input/combined_featurecounts.csv \
    data/input/samples.csv \
    data/processed/featurecounts.mat.RDS

# Build TPM table for other analyses
singularity exec -H $PWD:/home src/R.sif Rscript \
    src/buildTPMtable.R \
    data/input/combined_featurecounts.csv \
    data/input/samples.csv \
    data/processed/TPM.txt.gz

# Build DESeq Data Set (DDS) Object
singularity exec -H $PWD:/home src/R.sif Rscript \
    src/buildDDS.R \
    data/processed/featurecounts.mat.RDS \
    data/input/samples.csv \
    data/processed/DDS.RDS

# Run differential expression contrasts
# using src/runContrast.R

```


## Run DESeq
```
Rscript src/runDESeq.R \
    data/input/combined_featurecounts.csv \
    data/input/samples.csv \
    data/processed/featurecounts.mat.RDS
    data/processed/pseudohaploid_DEseq_dds.RDS
```

# Build heatmap
```
Rscript src/heatmap.R
```

# Run exploratory analyses
```
singularity run -H $PWD:/home src/R.sif
```# pseudo
