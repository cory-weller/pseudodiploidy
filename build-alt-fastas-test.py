#!/usr/bin/env python

import os
import sys
import subprocess
import gzip
from contextlib import ExitStack


chromosome = sys.argv[1]
start = sys.argv[2]
end = sys.argv[3]
chromosome = sys.argv[4]


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('vcf', type=str,
                    help='''Filename for vcf file. Can be gzipped ('.gz') or uncompressed ('.vcf')''')
parser.add_argument('reference', type=str,
                    help='''Filename for reference fasta file. Can be gzipped ('.gz') or uncompressed ('.fa', 'fasta')''')
parser.add_argument("--missing-case",
                    type=str,
                    nargs='?',
                    const=1,
                    default='lower',
                    help='''Accepts 'lower' and 'upper' as arguments.
                            If allele information is missing ('.') within vcf,
                            should the character be printed in <lower> (default) or
                            in <upper> case? See also: --missing-character''')
parser.add_argument("--missing-character",
                    type=str,
                    nargs='?',
                    const=1,
                    default='N',
                    help='''Accepts 'N', 'reference' or 'alternate' as arguments.
                            If allele information is missing ('.') within vcf,
                            should the character be printed as N (default) or
                            either reference or alternate alleles? If multiple
                            alternate alleles are available, one is chosen at random.
                            See also: --missing-case''')


args = parser.parse_args()

args.missing_character
args.missing_case
args.vcf
args.reference


class variantLibrary:
    help = 'Variant with reference and alternate allele(s)'
    def __init__(self, vcf, fasta) :
        self.db = {}
    def add(self, vcfLine):
        self.refStart = 
        self.refEnd = 
        self.refSeq = 
        self.altSeqs = []

strains = ['strain1', 'strain2']
from contextlib import ExitStack
variants = {
'strain1' : 'A',
'strain2' : 'G',
}

contig = 'chr1'
with ExitStack() as stack:
    for strain in strains:
        stack.enter_context(open(strain + '-' + contig + '.fasta', 'w'))
    for strain in strains]
    for i in outputConnections:
        i.write(variants[i.name])



vcfFileName = 'data/external/chromosome12.vcf.gz'
chromosomeName = 'chromosome12'

# Read through VCF and build list 'blocks' of genotypes
start, end, ref, alt1, alt2, alt3 ... etc

fasta = 'actgcaattgac'
positions =
51, 51 + len(ref)
<51+len(ref), next)
53
100
120
def printVariantFasta(vcf, fasta, strainsOfInterest):
    variantLibrary = {}
    variantLibrary[pos] = variant(arg1, arg2)
    genotypes = {}

    genotypes[posStart]

    prevPositionStart = 0
    prevPositionEnd = 0
    for line in vcf:
        splitLine = line.strip().split()[0:5]
        chromosome = splitLine[0]
        currentPositionStart = splitLine[1]
        assert currentPositionStart >= prevPositionEnd, 'current start position (%s) is less than previous ending position (%s)!' % (currentPositionStart, prevPositionEnd) 
        fullyReferenceRange = slice(prevPositionEnd, currentPositionStart)
        if currentPositionStart > prevPositionEnd:
            genotypes[posStart] = [fasta[fullyReferenceRange]]

        ref = splitLine[3]
        currentPositionEnd = currentPositionStart + len(ref)
        altGenotypes = splitline[4].split(',')
        strainGenotypes = splitline[9:]
        currentPositionEnd = currentPositionStart + len
        if currentPosition > prevPosition + len(splitLine[3]):
            genotypes[




# Reference allele = 

def buildGenotypeDictionary(vcfFileName):
    # Create dictionary for this chromosome
    genotypeDict = {}
    with gzip.open(args.vcf, 'r') as infile:
        for line in infile:
            try:                                                    # This try/except block allows for interfacing with gzip files;
                line = line.decode()                                # if the file is not gzipped, then it is opened like normal. If the
            except(UnicodeDecodeError, AttributeError):             # file is gzipped, then the line needs to be decoded from byte to string.
                pass                                                # A regular file throws UnicodeDecodeError; pass ignores this.
            chromosome, pos, _, ref, alt = line.strip().split()[0:5]
            strainGenotypes = line.strip().split()[9:]
            for strain, genotype in zip(strains, strainGenotypes):
                asldfkj
            if chromosome not in genotypeDict:
                genotypeDict[chromosome] = {}
            genotypes = [ref] + alt.split(',')      # index0 = ref, index1 = alt1, index2 = alt2, etc
            genotypeDict[chromosome][int(pos)] = genotypes
    return(genotypeDict)

genotypeDict = buildGenotypeDictionary(args.vcf)

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