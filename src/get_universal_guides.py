#!/usr/bin/env python


import sys
import argparse
import regex as re

def wrap_fasta(seq):
    '''wraps sequence every 80 characters with newlines'''
    return '\n'.join([seq[x:x+80] for x in range(0,len(seq),80)])

def import_vcf_header(filename):
    with open(filename, 'r', encoding='utf-8') as infile:
        header = infile.readlines()
        assert len(header) == 1, 'Expected header file with a single line. '\
                                  f'Offending file: {filename}'
        header = header[0].strip().split()
        return(header)



def subset_vcf(filename, chromosome, start, stop, strains='all'):
    with open(filename, 'r', encoding='utf-8') as infile:
        for line in infile:
            splitline = line.strip().split()
            if line.startswith('#'):
                if strains == 'all':
                    all_strains = header[9:]
                    strain_indices = [splitline.index(x) for x in all_strains]
                else:
                    strain_indices = [splitline.index(x) for x in strains]
                continue
            chr, pos = splitline[0:2]
            if pos < start:
                continue
            if pos >= start:
                return [splitline[x] for x in strain_indices]
            if pos > stop:
                break

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=str,
                    required = True,
                    help="""
                        Filename for fasta file. Can be gzipped ('.gz') 
                        or uncompressed ('.fa', 'fsa', 'fasta')
                    """
                    )
parser.add_argument('--start', type=int,
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

args = parser.parse_args()

if args.header_filename: # if header is provided
    print('TRUE')
else:
    print('FALSE')



header = import_vcf_header('./data/input/vcf-header.txt')

a = subprocess.run([
    "src/extract_fasta_range.py",
    "testref.fasta",
    '1',
    '10'
    ], capture_output=True)

a = ''.join(run([
    "src/extract_fasta_range.py",
    "testref.fasta",
    '1',
    '100'
    ], capture_output=True).stdout.decode().strip().split('\n')[1:])
a = ''.join(a)

a = subprocess.check_output([
    "src/extract_fasta_range.py",
    "testref.fasta",
    '1',
    '100'
    ], text=True)