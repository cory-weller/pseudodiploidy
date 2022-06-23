# Yeast Strain Barcoding

## Setup python environment
```
mamba create -n pseudo python=3.8
mamba activate pseudo
mamba install -q regex pyyaml pylint
mamba install -c conda-forge -c bioconda snakemake
```


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

## Check windows for strain uniqueness for first and second pass
```bash
# edit data/input/strains-to-start-with.txt
# which includes the starting set of possible strains to consider

chromosome=12
start=171320
end=172320
replicates=1
module load singularity

# Get Hygromycin-sensitive alpha-strains (pass 1)
awk 'NR > 1 && $3=="TRUE" && $5=="b" {print $1}' data/processed/strain-info.tsv > data/processed/hyg-sensitive-alpha.txt
src/query-region.sh data/processed/hyg-sensitive-alpha.txt ${chromosome} ${start} ${end} ${replicates} 
awk '$4=="b" {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | tr ',' '\n' > data/processed/hyg-sensitive-alpha-distinct.txt
rm reports/chr${chromosome}-${start}-${end}.tsv

# Get G418-sensitive a-strains (pass 1)
awk 'NR > 1 && $4=="TRUE" && $5=="a" {print $1}' data/processed/strain-info.tsv > data/processed/G418-sensitive-a.txt
src/query-region.sh data/processed/G418-sensitive-a.txt ${chromosome} ${start} ${end} ${replicates} 
awk '$4=="a" {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | tr ',' '\n' > data/processed/G418-sensitive-a-distinct.txt
rm reports/chr${chromosome}-${start}-${end}.tsv

# Do pass 2
src/get-second-pass-strains.py data/processed/hyg-sensitive-alpha.txt data/processed/hyg-sensitive-alpha-distinct.txt > data/processed/hyg-sensitive-alpha-pool2.txt
src/query-region.sh data/processed/hyg-sensitive-alpha-pool2.txt ${chromosome} ${start} ${end} ${replicates} 
awk '$4=="b" {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | tr ',' '\n' > data/processed/hyg-sensitive-alpha-pool2.txt
rm reports/chr${chromosome}-${start}-${end}.tsv

src/get-second-pass-strains.py data/processed/G418-sensitive-a.txt data/processed/G418-sensitive-a-distinct.txt > data/processed/G418-sensitive-a-pool2.txt
src/query-region.sh data/processed/G418-sensitive-a-pool2.txt ${chromosome} ${start} ${end} ${replicates} 
awk '$4=="a" {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | tr ',' '\n' > data/processed/G418-sensitive-a-pool2.txt
rm reports/chr${chromosome}-${start}-${end}.tsv


# extract strains identified in the first pass

# awk 'NR==2 {print $6}' reports/chr12-171320-172320.firstpass.tsv | sed 's/,/\n/g' | awk 'NR > 1' > data/input/first-pass-strains.txt
# awk 'NR==3 {print $6}' reports/chr12-171320-172320.firstpass.tsv | sed 's/,/\n/g' | awk 'NR > 1' >> data/input/first-pass-strains.txt

# python3 src/get-second-pass-strains.py

# src/query-region.sh data/input/second-pass-strains.txt ${chromosome} ${start} ${end} ${replicates} && \
# mv reports/chr${chromosome}-${start}-${end}.tsv reports/chr${chromosome}-${start}-${end}.secondpass.tsv
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
singularity exec --bind $(readlink ~/pseudodiploidy) src/R.sif Rscript \
    src/build-count-matrix.R \
    data/input/combined-featurecounts.csv \
    data/input/samples.csv \
    data/processed/featurecounts-matrix.RDS

# Build TPM table for other analyses
singularity exec --bind $(readlink ~/pseudodiploidy) src/R.sif Rscript \
    src/build-TPM-table.R \
    data/input/combined-featurecounts.csv \
    data/input/samples.csv \
    data/processed/TPM.txt.gz

# Build DESeq Data Set (DDS) Object
singularity exec --bind $(readlink ~/pseudodiploidy) src/R.sif Rscript \
    src/build-DDS.R \
    data/processed/featurecounts-matrix.RDS \
    data/input/samples.csv \
    data/processed/DDS.RDS

# Run differential expression contrasts
singularity exec --bind $(readlink ~/pseudodiploidy) src/R.sif Rscript \
    src/run-contrast.R
```

# Run exploratory analyses
```
singularity exec --bind $(readlink ~/pseudodiploidy) src/R.sif R
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
chr12:171320-172320

# look at guides across whole pool, first pass
cat data/processed/hyg-sensitive-alpha-distinct.txt data/processed/G418-sensitive-a-distinct.txt | sort -u > data/processed/pooled-distinct-strains.txt

singularity exec --bind $(readlink ~/pseudodiploidy) src/pseudodiploidy.sif Rscript src/get-variable-sites.R data/external/chromosome12.vcf.gz data/processed/pooled-distinct-strains.txt ${chromosome} ${start} ${end} > data/processed/pooled-distinct-strains-variable-sites.txt

singularity exec --bind $(readlink ~/pseudodiploidy) src/pseudodiploidy.sif Rscript src/get-usable-guides.R data/processed/chr12-guides.tsv data/processed/pooled-distinct-strains-variable-sites.txt $start $end


# same as above, but with second pool added
cat data/processed/G418-sensitive-a-pool2.txt data/processed/hyg-sensitive-alpha-pool2.txt data/processed/G418-sensitive-a-distinct.txt data/processed/hyg-sensitive-alpha-distinct.txt | sort -u > data/processed/double-pool-distinct.txt

singularity exec --bind $(readlink ~/pseudodiploidy) src/pseudodiploidy.sif Rscript src/get-variable-sites.R data/external/chromosome12.vcf.gz data/processed/double-pool-distinct.txt ${chromosome} ${start} ${end} > data/processed/pooled-distinct-strains-pool2-variable-sites.txt

singularity exec --bind $(readlink ~/pseudodiploidy) src/pseudodiploidy.sif Rscript src/get-usable-guides.R data/processed/chr12-guides.tsv data/processed/pooled-distinct-strains-pool2-variable-sites.txt $start $end

```R
library(data.table)
dat <- fread('data/external/chromosome12.vcf.gz')

start <- 171320
stop <- 172320

headerNames <- colnames(fread('data/input/vcf-header.txt', header=TRUE))
setnames(dat, headerNames)

usedStrains <- fread('data/processed/double-pool-distinct.txt', header=FALSE)$V1
usedStrains <- usedStrains[usedStrains!="S288C"]
setnames(dat, "#CHROM", "CHROM")

dat.sub <- dat[, c("CHROM", "POS", "REF", "ALT", usedStrains), with=FALSE][POS >= start & POS <= stop]
dat.sub[POS==171721]
```

# Looking at illumina seq data may 4 2022
Retrieve reads from nih box:
```
module load rclone
rclone copy --progress nihbox:/SBGE/data/simone/SMG_03142022/Fastq/8_S8_L001_R1_001.fastq.gz data/input/
rclone copy --progress nihbox:/SBGE/data/simone/SMG_03142022/Fastq/8_S8_L001_R2_001.fastq.gz data/input/

module load pear
pear -f data/input/8_S8_L001_R1_001.fastq.gz \
     -r data/input/8_S8_L001_R2_001.fastq.gz \
     -o data/processed/sample8

# remove unassembled or discarded reads, retaining 93.515 of original reads
rm data/processed/sample8.{discarded,unassembled.forward,unassembled.reverse}.fastq
gzip data/processed/sample8.assembled.fastq
```

# Making downsized VCF within window of interest
```R
library(data.table)
start <- 171320
end <- 172320
vcf <- fread('data/external/chromosome12.vcf.gz')
header <- colnames(fread('data/input/vcf-header.txt'))
setnames(dat, header)
strains <- fread('data/external/mating-types.tsv')[,Standardized_name]
wanted_cols <- c(header[1:9], strains)

vcf_window <- vcf[, ..wanted_cols][, !c("ID","QUAL","INFO","FILTER","FORMAT")][POS %between% c(start, end)]

fwrite(vcf_window, file='data/processed/window_genotypes.vcf', quote=F, row.names=F, col.names=T, sep='\t')

```