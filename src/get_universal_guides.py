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
            print(f'Importing header from file {header_filename}')
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

def gopen(filename):
    '''Generator that iterates through a file whether gzipped or not'''
    if filename.endswith('.gz'):
        open_fn = gzip.open
    else:
        open_fn = open
    with open_fn(filename, 'r') as infile:
        for line in infile:
            try:
                line = line.decode()
            except(UnicodeDecodeError, AttributeError):
                pass
            yield line


def subset_vcf(filename, chrom, start, stop, col_indices):
    '''loads in VCF for given chromosome, start/stop indices, for provided strains.
        Strains can be provided as list (loaded from file with --strains-file)'''
    for line in gopen(filename):
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
            yield [pos] + [splitline[x] for x in col_indices]

def as_string(*args, sep='\t'):
    '''converts multiple arguments to string, separated by <sep> (default of \t)'''
    return sep.join([str(x) for x in args])

def list_to_string(i, sep='\t'):
    '''converts list to string, separated by <sep> (default of \t)'''
    return sep.join([str(x) for x in i])

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

def get_fasta(filename, text_string='>'):
    '''returns the first entry in the fasta file that matches the provided text_string.
    If no text_string is provided, the first fasta entry will be used.'''
    found = False
    for line in gopen(filename):
        if found:
            if line.startswith('>') and text_string not in line:
                break
            yield line
        if text_string in line:
            found = True
            yield line




def format_fasta(fasta):
    '''formats the newline-separated fasta as a two-item list, [header,seq]'''
    fasta = list(fasta)
    if len(fasta) == 0:
        return None
    fasta = [x.strip() for x in fasta]
    return (fasta[0].lstrip('>'), ''.join(fasta[1:]))

def slice_fasta(fasta, start, stop):
    '''returns the fasta indices between start and stop (inclusive, 1-indexed)'''
    return fasta[(start-1) : stop]


def revcomp(seq):
    '''returns reverse complement of given sequence'''
    seq = seq.upper()
    pairs = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"
    }
    return ''.join([pairs[x] for x in seq[::-1]])

class guide:
    def __init__(self, match, strand, offset):
        self.offset = offset
        self.strand = strand
        if strand == '+':
            self.start = match.span()[0] + 1 + self.offset
            self.stop = (self.start + 19)
            self.seq = match.group(0)
            self.pam = self.seq[-3:]
            self.seq = self.seq[:20]
        elif strand == '-':
            self.start = match.span()[1] + self.offset
            self.stop = self.start - 19
            self.seq = revcomp(match.group(0))
            self.pam = self.seq[-3:]
            self.seq = self.seq[:20]

def findGuides(seq, forwardGuidePattern, reverseGuidePattern, offset):
    for forwardGuide in re.finditer(forwardGuidePattern, seq):
        yield guide(forwardGuide, '+', offset)
    for reverseGuide in re.finditer(reverseGuidePattern, seq):
        yield guide(reverseGuide, '-', offset)

parser = argparse.ArgumentParser()

parser.add_argument('--fasta', type=str,
                    required = True,
                    help="""
                        Filename for fasta file. Can be gzipped ('.gz') 
                        or uncompressed ('.fa', 'fsa', 'fasta')
                    """
                    )
parser.add_argument('--sequence-header', type=str,
                    required = False,
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
parser.add_argument('--out', type=str,
                    required = True,
                    help="""
                        filename (with or without absolute or relative path) for output files
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
parser.add_argument('--strains-list',
                    type=str,
                    nargs='?',
                    const=1,
                    help="""
                        Comma-separated string of strains to include
                    """
                    )
args = parser.parse_args()

# Import header from VCF itself, or provided filename)
header = import_vcf_header(args.vcf, args.header_filename)

# Import file name containing strains to use if provided
# or default to using all strains in VCF file if no file name provided
if args.strains_file:
    print(f'Importing header from {args.strains_file}')
    with open(args.strains_file, 'r', encoding='utf-8') as infile:
        strains = infile.readlines()
        strains = [x.strip() for x in strains]
elif args.strains_list:
    strains = args.strains_list.strip().split(',')
else:
    strains = header[9:]

print(f'VCF file: {args.vcf}')
print(f'Fasta file: {args.fasta}')
print('Desired Strains:')
for i in strains:
    print(i)

col_indices = [header.index(x) for x in strains]

print('Columns being retrieved from VCF:')
for i in col_indices:
    print(header[i])

fasta_header, fasta_seq = format_fasta(get_fasta(args.fasta, args.header_search))
fasta_seq = slice_fasta(fasta_seq, args.start, args.stop)

print(f'Header of FASTA entry selected: {fasta_header}')

#vcf = subset_vcf(filename, chrom, start, stop, strains):

vcf = subset_vcf(args.vcf, args.chrom, args.start, args.stop, col_indices)


print('Generating alternate fasta for selected strains...')
alt_fasta = list(fasta_seq)

with open(f'{args.out}.{args.start}-{args.stop}.variants', 'w', encoding='utf-8') as outfile:
    #outfile.write('POS\t' + list_to_string(strains) + '\n')
    for i in vcf:
        pos = str(i[0])
        genotypes = i[1:]
        variant = [strains[j] for j in range(1, len(genotypes)) if genotypes[j] == '1/1']
        if len(variant) == 0:
            continue
        outfile.write(pos +'\t' + str(len(variant)) + '\t' + '\t'.join(variant) + '\n')
        fixed_index = i[0]
        genotypes = set(i[1:])
        fixed_index = fixed_index - args.start
        if genotypes != {'0/0'}:
            alt_fasta[fixed_index] = 'N'
            #outfile.write(list_to_string(i) + '\n')

alt_fasta = ''.join(alt_fasta)
print('Done generating alternate fasta')

print('Finding for guide sequences compatible with selected strains...')
forwardGuidePattern = re.compile("[ACTG]{20}[ACTG]GG")
reverseGuidePattern = re.compile("CC[ACTG][ACTG]{20}")

guides = findGuides(alt_fasta, forwardGuidePattern, reverseGuidePattern, args.start)
n = 0
with open(f'{args.out}.guides.txt', 'w', encoding='utf-8') as outfile:
    outfile.write(as_string('guide', 'strand', 'PAM', 'start', 'end') + '\n')    
    for i in guides:
        n += 1
        outfile.write(as_string(i.seq, i.strand, i.pam, i.start, i.stop) +'\n')
print(f'Done, {n} guides found for {args.header_search} from {args.start}-{args.stop}')

# python src/get_universal_guides.py \
#     --out chr1 \
#     --vcf data/external/chromosome1.vcf.gz \
#     --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz \
#     --start 73986 \
#     --stop 74985 \
#     --header-filename data/input/vcf-header.txt \
#     --header-search "[chromosome=I]" \
#     --strains-file teststrains.txt

# python src/get_universal_guides.py \
#     --out chr2 \
#     --vcf data/external/chromosome2.vcf.gz \
#     --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz \
#     --start 137288 \
#     --stop 138287 \
#     --header-filename data/input/vcf-header.txt \
#     --header-search "[chromosome=II]" \
#     --strains-file teststrains.txt

# python src/get_universal_guides.py \
#     --out chr7 \
#     --vcf data/external/chromosome7.vcf.gz \
#     --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz \
#     --start 988749 \
#     --stop 989748 \
#     --header-filename data/input/vcf-header.txt \
#     --header-search "[chromosome=VII]" \
#     --strains-file teststrains.txt

# python src/get_universal_guides.py \
#     --out chr16 \
#     --vcf data/external/chromosome16.vcf.gz \
#     --fasta data/external/S288C_reference_sequence_R64-3-1_20210421.fsa.gz \
#     --start 618979 \
#     --stop 619978 \
#     --header-filename data/input/vcf-header.txt \
#     --header-search "[chromosome=XVI]" \
#     --strains-file teststrains.txt

