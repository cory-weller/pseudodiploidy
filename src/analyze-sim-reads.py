#!/usr/bin/env python


import itertools
import regex as re
import gzip
import bisect
import sys

#args = sys.argv
args = ['analyze-sim-reads.py', 'data/processed/sample8.assembled.fastq.gz', 'data/processed/sample8']
infile = args[1]
outStem = args[2]
goodSeqFile = '%s.pass.txt' % outStem
badSeqFile = '%s.fail.txt' % outStem


constantL = [   'GCGAAATGGGTAATCTCTTTACAGTTAACAGATTCTCTGGTGATTGGCTT', 
                'ATATCCGGTCTCGATTTCTAATAACGGAGATTTTGTTGGCCAAGAGACA',
                'TGACAGAATATGCCAAAGAACCCATAAATAAATATGATATAAGAGC'
            ]

constantR = [   'CTCAAGGTAGCAAATATACTCTAATTGTAGCACGTCTACTATGTATTT',
                'AAGGGCAGCACCAAAGTTAAGTAAAAAGAGACTTTCAAATACATTGGAT',
                'CACTGGGCCGGCGTTGGTCAGAGGTGTGGATAAACCAATGAAAAGACCTGT'
            ]

constantTolerance = 0

patternL = re.compile(r"""\L<string>{s<=%s}""" % (constantTolerance), string=constantL)
patternR = re.compile(r"""\L<string>{s<=%s}""" % (constantTolerance), string=constantR)


def revcomp(seq):
    pairs = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C",
        "N" : "N"
    }
    return ''.join([pairs[x] for x in seq.upper()[::-1]])

def splitSeq(myseq):
    # Trim left constant sequence
    matchL = re.search(patternL, myseq).captures()[0]
    try: beforeMatchL, afterMatchL = myseq.split(matchL)
    except ValueError: return(None)
    # Trim right plasmid sequence from everything after the left match
    matchR = re.search(patternR, afterMatchL).captures()[0]
    try: beforeMatchR, afterMatchR = afterMatchL.split(matchR)
    except ValueError: return(None)
    return([beforeMatchL, matchL, beforeMatchR, matchR, afterMatchR])

print("Importing gzipped fastq sequences from %s" % infile)


with gzip.open(infile, 'r') as f, open(goodSeqFile, 'w') as goodOut, open(badSeqFile, 'w') as badOut:
    for lines in itertools.zip_longest(*[f]*4):      
        try:
            rawSeq = lines[1].strip().decode()
            result = splitSeq(rawSeq)
            goodOut.write('\t'.join(result) + '\n')
        except(UnicodeDecodeError, AttributeError):
            badOut.write(rawSeq + '\n')


library(data.table)
dat <- fread('data/processed/sample8.pass.txt', header=F)
dat[, 'V1' := NULL]
dat[, 'V5' := NULL]

dat[, .N, by=list(V2, V3)]