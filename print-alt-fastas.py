#!/usr/bin/env python

import os
import sys
import subprocess
import gzip

chromosome = sys.argv[1]
start = sys.argv[2]
end = sys.argv[3]
chromosome = sys.argv[4]

vcfFileName = 'data/external/chromosome12.vcf.gz'
chromosomeName = 'chromosome12'

def buildGenotypeDictionary(vcfFileName):
    # Create dictionary for this chromosome
    genotypeDict = {}
    with gzip.open(vcfFileName, 'r') as infile:
        for line in infile:
            try: line = line.decode()
            except(UnicodeDecodeError, AttributeError): pass
            chromosome, pos, _, ref, alt = line.strip().split()[0:5]
            strainGenotypes = line.strip().split()[9:]
            for strain, genotype in zip(strains, strainGenotypes):
                asldfkj
            if chromosome not in genotypeDict:
                genotypeDict[chromosome] = {}
            genotypes = [ref] + alt.split(',')      # index0 = ref, index1 = alt1, index2 = alt2, etc
            genotypeDict[chromosome][int(pos)] = genotypes
    return(genotypeDict)

genotypeDict = buildGenotypeDictionary(vcfFileName)
positions = list(genotypeDict[chromosomeName].keys())
positions.sort(reverse = True)

# fasta = ''
# import VCF as list
# sort inverse by position (descending, high to low)
# for item in list, 
# read in whole thing
# reverse order
# 

vcf[position][individual]
def import_fasta(filename):
    with gzip.open(filename, 'r') as infile:
        text = infile.read().splitlines()
        text = list(''.join(text[1:]))
    return(text)

def personalize_fasta(site_genotypes, personal_genotypes_file, chromosome, ref_sequence):
    haplotype1 = ref_sequence[:]
    haplotype2 = ref_sequence[:]
    with open(personal_genotypes_file, 'r') as genotype_file:
        for line in genotype_file:
            chrom, pos, hap1, hap2 = line.split()
            if chrom == chromosome:
                pos = int(pos)
                if hap1 == "1/1":
                    haplotype1[pos-1] = site_genotypes[pos][1]
                elif hap1 == "./.":
                    haplotype1[pos-1] = "N"
                if hap2 == "1/1":
                    haplotype2[pos-1] = site_genotypes[pos][1]
                elif hap2 == "./.":
                    haplotype2[pos-1] = "N"
    return(haplotype1, haplotype2)

def write_fasta(ind_n, chromosome, haplotype1, haplotype2, wrap_length):
        with open(str(ind_n) + "." + chromosome + ".fasta", 'w') as outfile:
          outfile.write(">" + str(ind_n) + "_" + chromosome + "_haplotype1" "\n")
          outfile.write('\n'.join([''.join(haplotype1[i:i+wrap_length]) for i in range(0,len(haplotype1),wrap_length)]))
          outfile.write('\n')
          outfile.write(">" + str(ind_n) + "_" + chromosome + "_haplotype2" "\n")
          outfile.write('\n'.join([''.join(haplotype2[i:i+wrap_length]) for i in range(0,len(haplotype2),wrap_length)]))

def get_diploidGenomeSize(fastaFilename):
    cmd = """grep "^[^>]" %s | wc -c""" % (fastaFilename)
    output = subprocess.check_output(cmd, shell=True).rstrip()
    return(int(output))

def run_wgsim(fastaFilename, readLength, coverage):
    diploidGenomeSize = get_diploidGenomeSize(fastaFilename)
    nReads = round( (float(diploidGenomeSize) * coverage) / (4*readLength) )
    ind_n, chromosome = fastaFilename.split(".")[0:2]
    cmd = """../../bin/wgsim-master/wgsim \
    -1 %s \
    -2 %s \
    -N %s \
    -e 0.001 \
    -r 0 \
    -R 0 \
    %s.%s.fasta \
    %s.%s.F.fq \
    %s.%s.R.fq && \
    rm %s.%s.fasta && \
    bzip2 %s.%s.F.fq && \
    bzip2 %s.%s.R.fq""" % (readLength, readLength, nReads, ind_n, chromosome, ind_n, chromosome, ind_n, chromosome, ind_n, chromosome, ind_n, chromosome, ind_n, chromosome)
    subprocess.call(cmd, shell=True)

args = sys.argv[1:]
project_directory = os.path.dirname(os.getcwd())

population = args[0]
chromosome = args[1]
first_ind = int(args[2])
last_ind = int(args[3])

site_genotypes = build_ref_alt_dict(directory = project_directory + "/input_data/", chromosome=chromosome)

ref_sequence = import_fasta(filename = project_directory + "/input_data/" + chromosome + ".fa")

os.chdir(population)

# Iterate through individuals
for ind_n in range(first_ind, last_ind+1, 1):
    fastaFilename = str(ind_n) + "." + chromosome + ".fasta"
    if os.path.exists(str(ind_n) + "." + chromosome + ".R.fq.bz2") == False:
        haplotype1, haplotype2 = personalize_fasta( site_genotypes=site_genotypes,
                                                    personal_genotypes_file=str(ind_n) + ".geno",
                                                    chromosome=chromosome,
                                                    ref_sequence=ref_sequence)
        write_fasta(ind_n, chromosome, haplotype1, haplotype2, wrap_length=80)
        run_wgsim(fastaFilename, readLength=100.0, coverage=0.05)