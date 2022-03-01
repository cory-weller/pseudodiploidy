# Yeast Pseudodiploidy

## Singularity
### Install Singularity
```
bash src/install-singularity.sh
```

### Build Singularity image
```
singularity build --fakeroot src/pseudodiploidy.sif src/pseudodiploidy.def
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
singularity exec --bind ${PWD} src/R.sif Rscript \
    src/run-contrast.R
```

# Run exploratory analyses
```
singularity exec --bind ${PWD} src/R.sif R
```

# 1011-barcodes


Retrieved table of genes and gene predictions for S. cerevisiae from UCSC Table Browser

Search settings for gene coordinates:

| Field     | Value                                 |
| -----     | -----                                 |
| clade     | other                                 |
| genome    | S. cerevisiae                         |
| assembly  | Apr. 2011 (SacCer_Apr2011/sacCer3)    |
| group     | Genes and Gene Predictions            |
| track     | SGD Genes                             |
| table     | sgdGene                               |
| region    | genome                                |


Search settings for expression coordinates:

| Field     | Value                                 |
| -----     | -----                                 |
| clade    | other                                  |
| genome    | S. cerevisiae                         |
| assembly  | Apr. 2011 (SacCer_Apr2011/sacCer3)    |
| group     | Expression and Regulation             |
| track     | Regulatory Code                       |
| table     | transRegCode                          |
| region    | genome                                |

# retrieve bgzipped and tabix-index for S288C vcf files
```
remotedata='nihbox'
mkdir -p data/genomes/
rclone copy ${remotedata}:/S288C/vcf/ data/genomes/
# If working with raw vcf files, generate bgzipped files
# and tabix indices via:
# (cd data/genomes && for file in *.vcf; do singularity exec ../src/singularity.sif bgzip $file; done)
# (cd data/genomes && for file in *.vcf.gz; do singularity exec ../src/singularity.sif tabix -p vcf $file; done)

```

# extract regions from vcf
singularity exec src/singularity.sif tabix genomes/chromosome1.vcf.gz chromosome1:1-100

# test region

singularity exec src/singularity.sif Rscript stuff.R 1 500 1000