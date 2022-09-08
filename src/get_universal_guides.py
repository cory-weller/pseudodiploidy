#!/usr/bin/env python


import sys
import argparse
import subprocess
import gzip
import regex as re

def wrap_fasta(seq):
    '''wraps sequence every 80 characters with newlines'''
    return '\n'.join([seq[x:x+80] for x in range(0,len(seq),80)])

def import_vcf_header(vcf_filename, header_filename=None):
    '''imports (separate) file containing VCF header names'''
    vcf_header = None
    external_header = None
    if header_filename:
        with open(header_filename, 'r', encoding='utf-8') as infile:
            external_header = infile.read().strip().split()
    with gzip.open(vcf_filename) as infile:
        firstline = infile.readline().decode()
        if firstline.startswith('#'):
            vcf_header = firstline.strip().split()
        if external_header and vcf_header:
            assert vcf_header == external_header, 'Provided header does not match VCF header!'
            return external_header
    if vcf_header and not external_header:
        return vcf_header
    if external_header and not vcf_header:
        return external_header



def subset_vcf(filename, chrom, start, stop, col_indices):
    '''loads in VCF for given chromosome, start/stop indices, for provided strains.
        Strains can be provided as list (loaded from file with --strains-file'''
    with gzip.open(filename) as infile:
        for line in infile:
            line = line.decode()
            splitline = line.strip().split()
            if line.startswith('#') and header is None:
                continue
            chrom, pos = splitline[0:2]
            pos = int(pos)
            if pos < start:
                continue
            if pos > stop:
                break
            if pos >= start:
                yield [splitline[x] for x in col_indices]


def extract_fasta(fasta_filename, start, stop, strip_header=False, header_search=None):
    '''returns given start/stop range for provided fasta file'''
    args = [ 'src/extract_fasta_range.py',
                                    f'{fasta_filename}',
                                    f'{start}',
                                    f'{stop}']
    if header_search:
        args.append('--header-search')
        args.append(header_search)
    out = subprocess.check_output(args, text=True)
    if strip_header:
        out = ''.join(out.strip().split('\n')[1:])
    return out


parser = argparse.ArgumentParser()

parser.add_argument('--fasta', type=str,
                    required = True,
                    help="""
                        Filename for fasta file. Can be gzipped ('.gz') 
                        or uncompressed ('.fa', 'fsa', 'fasta')
                    """
                    )
parser.add_argument('--vcf', type=str,
                    required = True,
                    help="""
                        Filename for VCF file
                    """
                    )
parser.add_argument('--chrom', type=str,
                    required = False,
                    help="""
                        Chromosome for VCF subsetting. Currently not yet implemented / used
                    """
                    )
parser.add_argument('--start', type=int,
                    required = True,
                    help="""
                        1-indexed start position to begin with (inclusive)
                    """
                    )
parser.add_argument('--stop', type=int,
                    required = True,
                    help="""
                        1-indexed stop position to end at (inclusive)
                    """
                    )
parser.add_argument('--header-filename',
                    type=str,
                    nargs='?',
                    const=1,
                    default=None,
                    help="""
                        Name of single-line external file that defines VCF column headers
                    """
                    )
parser.add_argument('--header-search',
                    type=str,
                    nargs='?',
                    const=1,
                    default=None,
                    help="""
                        Pattern to search for in headers in order to determine which
                        sequence to use for extracting the range
                    """
                    )
parser.add_argument('--header-output',
                    type=str,
                    nargs='?',
                    const=1,
                    help="""
                        Name to include in output header, followed by <start>-<stop>
                    """
                    )
parser.add_argument('--strains-file',
                    type=str,
                    nargs='?',
                    const=1,
                    help="""
                        Strains (columns) to load in VCF file, one strain per line
                    """
                    )

args = parser.parse_args()

header = import_vcf_header(args.vcf, args.header_filename)

if args.strains_file:
    with open(args.strains_file, 'r', encoding='utf-8') as infile:
        strains = infile.readlines()
        strains = [x.strip() for x in strains]
else:
    strains = header[9:]

col_indices = [strains.index(x) for x in strains]


fasta = extract_fasta(args.fasta, args.start, args.stop, strip_header=True)
print(fasta)

#vcf = subset_vcf(filename, chrom, start, stop, strains):

vcf = subset_vcf(args.vcf, args.chrom, args.start, args.stop, col_indices)

#for i in vcf:
#    print(i)

# python src/get_universal_guides.py \
#     --vcf data/external/chromosome1.vcf.gz \
#     --fasta testref.fasta \
#     --start 10000 \
#     --stop 11000 \
#     --header-filename data/input/vcf-header.txt \
#     --strains-file teststrains.txt

# python src/get_universal_guides.py \
#     --vcf data/external/chromosome1.vcf.gz \
#     --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa \
#     --start 1 \
#     --stop 11000 \
#     --header-filename data/input/vcf-header.txt \
#     --header-search "chromosome=I" \
#     --strains-file teststrains.txt
