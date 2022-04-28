# Yeast Strain Barcoding

## Activate OR install Singularity
```
# install on your system if you have sudo
bash src/install-singularity.sh

# activate module, e.g. on HPC
module load singularity
```

## Retrieve data
```
bash src/get-data.sh
```

## Check windows for strain uniqueness
```
# run as
# src/query-region.sh <strains-to-use> <chromosome> <startPos> <endPos> <nReplicates>
# e.g.:
src/query-region.sh data/input/strains-to-start-with.txt 12 171320 172320 10

# results written to reports/ as chr#-start-stop.tsv
```

## Print regions of interest
```
# print to terminal
bash src/print-region.sh 2 80000 81000

# print to file
bash src/print-region-to-file.sh 2 80000 81000 

```


# Yeast Pseudodiploidy
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


# Finding putative guides
Retrieve S288C Reference genome
```
rclone copy nihbox:/cloud/S288C/S288C_reference_sequence_R64-3-1_20210421.fsa.gz data/external/ && \
gunzip data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz
```

Extract just chromosome 12
```
grep -n 'chromosome=XII' data/external/S288C_reference_sequence_R64-3-1_20210421.fsa
```
Which returns output:
```
120768:>ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]
138739:>ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]
```

So we want lines 120768:138738 (17971 lines):
```
head -n 138738 data/external/S288C_reference_sequence_R64-3-1_20210421.fsa | tail -n 17971 > data/external/S288C-chr12.fasta
```

Extract variable sites:
```
zcat data/external/chromosome12.vcf.gz | awk '{print $2}' > data/processed/S288C-chr12-variable-sites.txt
```
Find all guides
```
python3 src/find-guides.py > data/processed/chr12-guides.tsv
```
chr12:171320-172320

singularity exec ../pseudodiploidy/src/pseudodiploidy.sif Rscript src/get-variable-sites.R data/external/chromosome12.vcf.gz data/input/type-a-strains.txt 12 171320 172320 > varSites-a.txtR

vcf         data/external/chromosome12.vcf.gz
strains     data/input/type-a-strains.txt
strains     data/input/type-b-strains.txt
chromosome  12