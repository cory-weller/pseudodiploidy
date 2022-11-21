
```

src/query-region.sh data/input/haploid-strains.txt 16 618979 619978 1

48 strains first pass chr16:618979-619978

a   AFH,AHG,AGR,AFQ,AFD,AHQ,AFF,AHM,AAG,AEG,AHL,AIE,AGS,AFR,AGP,AHI,AHV,AFN,AGT
b   AEM,ADC,ACD,ABB,ABT,AEF,ADI,ABH,ADF,ADE,ACH,ABI,ACV,ABL,ABM,ADA,ACQ,ADD,ABC,ACR,AGC,ADG,ABF,ADR,ABG,ADV,ABP

chromosome=16; start=618979; end=619978; replicates=1
awk 'NR==2 {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | sed 's/,/\n/g' | awk 'NR > 1' > data/input/${chromosome}-${start}-${end}-firstpass.txt
awk 'NR==3 {print $6}' reports/chr${chromosome}-${start}-${end}.tsv | sed 's/,/\n/g' | awk 'NR > 1' >> data/input/${chromosome}-${start}-${end}-firstpass.txt


mkdir -p reports/chosen/
cat reports/chr${chromosome}-${start}-${end}.tsv > reports/chosen/chr${chromosome}-${start}-${end}.first.tsv

python3 src/get-second-pass-strains.py data/input/haploid-strains.txt data/input/${chromosome}-${start}-${end}-firstpass.txt > data/input/${chromosome}-${start}-${end}-secondpass-start.txt
src/query-region.sh data/input/${chromosome}-${start}-${end}-secondpass-start.txt ${chromosome} ${start} ${end} ${replicates}
mv reports/chr${chromosome}-${start}-${end}.tsv reports/chosen/chr${chromosome}-${start}-${end}.second.tsv

all chosen strains
a1  AFH,AHG,AGR,AFQ,AFD,AHQ,AFF,AHM,AAG,AEG,AHL,AIE,AGS,AFR,AGP,AHI,AHV,AFN,AGT
b1  AEM,ADC,ACD,ABB,ABT,AEF,ADI,ABH,ADF,ADE,ACH,ABI,ACV,ABL,ABM,ADA,ACQ,ADD,ABC,ACR,AGC,ADG,ABF,ADR,ABG,ADV,ABP
a2  AIA,AHA,AEH,ADN,AHH
b2  ABS,ACT,AAT,ACM,AIB,AFM,ACL,ACG,ACC
```
```bash 
usedstrains="AGL,ABE,ABR,ABQ,AHG,AHD,AGF,AHE,ACP,AHR,AGH,ACF,AAI,AHT,AAA,ADR,ADC,AGV,AGS,ADB,ABT,ADI,ADD,ACR,AHI,ABS,BIQ,CEA,CEB,AHN,ACT,ACV,ABP,AAK,AET,AEK,AHK,AES,AFI,AIA,AAE,AIB,AIC,AVB,ABM,AHM,AHL,AIE,AGN,AHF,ACS,ACG,ACH,AFR,AEG,AGM,AEH,AFD,AHC,APG,BVE,AHB,AFL,AAB,AFB,AAR,AEA,BRM,AEQ,AAQ,AIG,AFQ,AER,AFF,AAM,AHV,AEM,AGK,BCN,ADE,AAH,AGC,ACN,CBK,ADA,ACA,ACB,BRL,ADT,ABH,ABI,AAV,ABK,ABA,AAS,AAT,ABG,AGP,AAG,AHA,AGT,ABF,AEI,AFG,ABC,CDK,CEM,ACI,ACL,AHS,ADV,ADF,ACQ,ABL,ACC,ACD,ADG,AHH,CDA,AGR,AEF,CLQ,AFN,ADQ,AFM,ADN,AHQ,AFH,CEV,ACM,ABB,ACK,ADH"

IFS=',' read -r -a strainarray <<< "$usedstrains"



start=677227
end=677246
let start=${start}-50
let end=${end}+53
chr="XVI"
windowFasta="data/input/S288C_chr16_${start}-${end}.fa"

python src/extract_fasta_range.py \
    data/external/S288C_reference_sequence_R64-3-1_20210421.fsa \
    ${start} \
    ${end} \
    --header-search ${chr} > ${windowFasta}

# Generate BLAST database for comparing to individual strain queries
module load blast
makeblastdb -dbtype nucl -in ${windowFasta}

# iterate through strains perform BLAST search
mkdir -p tmp/
for strain in ${strainarray[@]}; do
    # unzip strain-specific fasta
    unzip -p data/external/1011-assemblies.zip ${strain}.fasta > tmp/${strain}.fasta
    blastn -db ${windowFasta} -outfmt 6 -query tmp/${strain}.fasta  > tmp/${strain}.blast
    while read query db pct length i j qstart qend dbstart dbend m; do
        if [[ $dbstart > $dbend ]]; then
            startpos=$qend
            endpos=$qstart
        else
            startpos=$qstart
            endpos=$qend
        fi
        python src/extract_fasta_range.py \
            tmp/${strain}.fasta \
            ${startpos} \
            ${endpos} \
            --header-search ${query} >> chr16-windows.fasta
    done < tmp/${strain}.blast
    rm tmp/${strain}.fasta
done &

python src/print-clustal.py chr16-window.msa > chr16-window.msa.txt

```






```
python src/get_universal_guides.py \
    --out chr16 \
    --vcf data/external/chromosome16.vcf.gz \
    --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz \
    --start 618979 \
    --stop 619978 \
    --header-filename data/input/vcf-header.txt \
    --header-search "[chromosome=XVI]" \
    --strains-list AFH,AHG,AGR,AFQ,AFD,AHQ,AFF,AHM,AAG,AEG,AHL,AIE,AGS,AFR,AGP,AHI,AHV,AFN,AGT,AEM,ADC,ACD,ABB,ABT,AEF,ADI,ABH,ADF,ADE,ACH,ABI,ACV,ABL,ABM,ADA,ACQ,ADD,ABC,ACR,AGC,ADG,ABF,ADR,ABG,ADV,ABP
```