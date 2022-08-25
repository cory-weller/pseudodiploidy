#!/usr/bin/env python

import sys
import gzip

chromosome = sys.argv[1]
start = sys.argv[2]
end = sys.argv[3]
chromosome = sys.argv[4]

vcf_filename = f"data/external/chromosome{chromosome}.vcf.gz"


def build_genotype_dictionary(vcf_filename):
    # Create dictionary for this chromosome
    #genotype_dict = {}
    pos_array = []
    ref_array = []
    alt_array = []
    with open(vcf_filename, 'r', encoding='utf-8') as infile:
        for line in infile:
            if line.startswith('#CHROM'):
                continue
            chromosome, pos, ref, alt = line.strip().split()[0:4]
            pos_array.append(int(pos))
            ref_array.append(ref)
            alt_array.append(alt)
    length_array = [len(x) for x in ref_array]
    next_pos_array = pos_array[1:] + [pos_array[-1]+length_array[-1]]
    ref_range_array = [range((pos + L),nextPos) for pos, L, nextPos in zip(pos_array, length_array,next_pos_array)]
    for i in zip(pos_array, ref_array, alt_array, length_array, next_pos_array, ref_range_array):
        print(i)
    #return(genotype_dict)

build_genotype_dictionary('data/processed/window_genotypes.vcf')


genotype_dict = buildgenotype_dictionary(vcf_filename)
positions = list(genotype_dict[chromosome_name].keys())
positions.sort(reverse = True)

# fasta = ''
# import VCF as list
# sort inverse by position (descending, high to low)
# for item in list, 
# read in whole thing
# reverse order
# 

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
