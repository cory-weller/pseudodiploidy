# README

This directory tests the validation of MAT-a and MAT-alpha classification, tested on the
parental lines used in the Josh Bloom 16 parental strains with known mating types.


```bash
PRJ='PRJNA549760'
esearch -db sra -query ${PRJ} | efetch -format runinfo > ${PRJ}_info.csv

download_sra() {
    local srr=${1}
    wget -O ${srr}.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/${srr}/${srr}
}

fasterq-dump ${srr}
```

```bash
< ${PRJ}_info.csv tr ',' '\t'   | \
    cut -f 1,30                 | \
    grep bam                    | \
    sed 's/\.bam$//g' \
> ${PRJ}_founders.txt
```

fasterq-dump SRR9330834.sra

```bash
module load sratoolkit
while read SRA STRAIN; do
    echo ${STRAIN} ${SRA}
done < ${PRJ}_founders.txt

download_sra SRR9330834
fasterq-dump SRR9330834.sra
