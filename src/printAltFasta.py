#!/usr/bin/env python

import sys
import gzip

# chromosome = sys.argv[1]
# start = sys.argv[2]
# end = sys.argv[3]
# chromosome = sys.argv[4]

# vcf_filename = sys.argv[1]
# fasta_filename = sys.argv[2]
# start = int(sys.argv[3])
# end = int(sys.argv[4])
# strain = sys.argv[5]

vcf_filename = 'data/external/chromosome12.vcf.gz'
header_filename = 'data/input/vcf-header.txt'
start = 171320
end = 172320
fasta_filename = 'data/external/S288C-chr12.fasta'
strain = sys.argv[1]
#strain = 'AAI'
wrap_length = 80

with open(fasta_filename, 'r', encoding='utf-8') as infile:
    seq = list('N' + ''.join([x.strip() for x in infile.readlines()[1:]]))
    seq = ['N'] + seq[(start) : (end+1)]

output = []
variants = 0

seq2 =''.join(seq[1:])
print('\n'.join([seq2[i:i+wrap_length] for i in range(0,len(seq2),wrap_length)]))
last_pos = 1
if strain != 'REF':
    with open(header_filename, 'r') as infile:
        for line in infile:
            if line.startswith('#CHROM'):
                splitline = line.strip().split()
                strain_idx = splitline.index(strain)

overlap = FALSE

with gzip.open(vcf_filename, 'r') as infile:
    for line in infile:
        constant_region = ''
        splitline = line.decode().strip().split()
        chromosome = splitline[0]
        raw_pos = int(splitline[1])
        if raw_pos < start:
            continue
        if raw_pos > end:
            break
        variants += 1
        pos = raw_pos - start + 1
        ref = splitline[3]
        alts = splitline[4].split(',')
        vcf_genotypes = [ref] + alts
        # NOTE: Takes only the second allele if heterozygous
        if strain == 'REF':
            strain_allele = 0
        else:
            strain_allele = splitline[strain_idx].split('/')[1]
            if strain_allele == '.':
                strain_allele = 0
            else:
                strain_allele = int(strain_allele)
        constant_region = seq[last_pos : pos]
        genotype = list(vcf_genotypes[strain_allele])
        #print(constant_region, ' ', genotype)
        if len(constant_region) > 0:
            output += constant_region
        output += genotype
        last_pos = pos + len(ref)
        #print(len(vcf_genotypes[0]))
        print(pos, len(output), vcf_genotypes, strain_allele)
        #print(''.join(output))
        #print(''.join(pre) + '\t' + ''.join(genotype) + '\t' +''.join(post))

print(variants)
seq = ''.join(output)
# offset works before genotype, need to update after
wrap_length = 80
print(f">{strain}_{chromosome}_{start}-{end}")
print('\n'.join([seq[i:i+wrap_length] for i in range(0,len(seq),wrap_length)]))
