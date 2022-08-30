#!/usr/bin/env python
'''finds a header within a given fasta file, then prints the desired nucleotide ranges to STDOUT'''

import sys
import argparse
import regex as re

def wrap_fasta(seq):
    '''wraps sequence every 80 characters with newlines'''
    return '\n'.join([seq[x:x+80] for x in range(0,len(seq),80)])

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str,
                    help="""
                        Filename for fasta file. Can be gzipped ('.gz') 
                        or uncompressed ('.fa', 'fsa', 'fasta')
                    """
                    )
parser.add_argument('start', type=int,
                    help="""
                        1-indexed start position to begin with (inclusive)
                    """
                    )
parser.add_argument('stop', type=int,
                    help="""
                        1-indexed stop position to end at (inclusive)
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

if args.header_search is None:
    header_pattern = '>'
else:
    header_pattern = args.header_search + r"[^IVX]+"
    header_pattern = re.compile(header_pattern)

# format header
fasta_filename = args.fasta.split("/")[-1]
out_header = re.split('.fasta|.fa|.fsa', fasta_filename)[0]
if args.header_search:
    out_header += '_' + args.header_search
out_header += ':'+ str(args.start) + '-' + str(args.stop)

out_seq = ''
store_seq = False
found_headers = []

with open(args.fasta, 'r', encoding='utf-8') as infile:
    for line in infile:
        if line.startswith('>'):
            if bool(re.search(header_pattern, line)):
                store_seq = True
                found_headers.append(line.strip())
            else:
                store_seq = False
            if len(found_headers) > 1:
                print("Found multiple matching headers:")
                print("\n".join(found_headers))
                print("Exiting")
                sys.exit()
            continue
        if store_seq is True:
            out_seq += line.strip()


if len(found_headers) == 0:
    print("no matches found")
    sys.exit()

if len(out_seq) == 0:
    print("seq is of length 0")
    sys.exit()

wanted_region = out_seq[(args.start-1) : args.stop]


if len(wanted_region) == 0:
    print("extracted region is of length 0")
    sys.exit()


if args.header_output:
    print('>' + args.header_output)
else:
    print('>' + out_header)
print(wrap_fasta(wanted_region))
