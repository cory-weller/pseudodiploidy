# Readme for finding new potential windows

for strains 
MSY32 and MSY33

 
# While in folder with pacbio_fasta.zip:
```
function splitFasta() {
    inFasta=${1}
    awk -F "(^>|\t| )" '{if($0 ~ /^>/) {s=$2".fasta";  print ">"$2 > s} else print > s}' ${inFasta}
}

export -f splitFasta
```


unzip -p pacbio_fasta.zip MSY32.fasta > MSY32.fasta
unzip -p pacbio_fasta.zip MSY33.fasta > MSY33.fasta

```targets.txt
I   71786   76152
XVI 615379  621258
VII 988049  993521
II  136688  140260
```

```bash
while read chromosome start stop; do
    python ~/pacbio-yeast-genomes/extract_fasta_range.py \
    ../data/external/S288C_reference_sequence_R64-3-1_20210421.fsa \
    ${start} \
    ${stop} \
    --header-search chromosome=${chromosome}
done < targets.txt > targets.fasta
```



module load blast
makeblastdb -dbtype nucl -in MSY32.fasta
blastn -db MSY32.fasta -outfmt 6 -query targets.fasta | cut -f 1,2,9,10 > MSY32.hits
while read ref query start stop; do
    chr=$(echo $query | cut -d "_" -f 2)
    python ~/pacbio-yeast-genomes/extract_fasta_range.py \
        MSY32.fasta \
        ${start} \
        ${stop} \
        --header-search ${chr}
done < MSY32.hits > MSY32.windows.fasta

makeblastdb -dbtype nucl -in MSY33.fasta
blastn -db MSY33.fasta -outfmt 6 -query targets.fasta | cut -f 1,2,9,10 > MSY33.hits
while read ref query start stop; do
    chr=$(echo $query | cut -d "_" -f 2)
    python ~/pacbio-yeast-genomes/extract_fasta_range.py \
        MSY33.fasta \
        ${start} \
        ${stop} \
        --header-search ${chr}
done < MSY33.hits > MSY33.windows.fasta

# output:
# chrXII  MSY29_chrXII    99.778  4501    4       1       1       4501    161646  166140  0.0     8251
# chrXII  MSY32_chrXII    99.067  4501    27      3       1       4501    160290  164775  0.0     8065
# chrXII  MSY33_chrXII    99.401  4509    5       2       1       4501    156082  160576  0.0     8157
# chrXII  MSY39_chrXII    99.134  4501    16      3       1       4501    165006  169483  0.0     8074


```python
#!/usr/bin/env python
'''extractRegion.py'''

import sys

fasta_name = sys.argv[1]
strain_name = fasta_name.split(".fasta")[0]
start_pos = int(sys.argv[2])
end_pos = int(sys.argv[3])
out_header = f">{strain_name}_{start_pos}-{end_pos}"

def wrapFasta(seq):
    return '\n'.join([seq[x:x+80] for x in range(0,len(seq),80)])

with open(fasta_name, 'r', encoding='utf-8') as infile:
    seq = infile.readlines()[1:]

seq = [x.strip() for x in seq]
# pad at start to account for python 0-indexing
seq = 'N' + ''.join(seq)
wanted_region = seq[start_pos : end_pos+1]

print(out_header)
print(wrapFasta(wanted_region))
```

python extractRegion.py MSY29_chrXII.fasta 161646 166140 > MSY29_chrXII_window.fasta
python extractRegion.py MSY32_chrXII.fasta 160290 164775 > MSY32_chrXII_window.fasta
python extractRegion.py MSY33_chrXII.fasta 156082 160576 > MSY33_chrXII_window.fasta
python extractRegion.py MSY39_chrXII.fasta 165006 169483 > MSY39_chrXII_window.fasta
