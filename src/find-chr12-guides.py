#!/usr/bin/env python3

import re

def revcomp(seq):
    pairs = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"
    }
    return ''.join([pairs[x] for x in seq.upper()[::-1]])

with open('data/external/S288C-chr12.fasta', 'r') as infile:
    fasta = ''.join([x.strip() for x in infile.readlines()[1:]])

fasta = list(fasta)

with open('data/processed/S288C-chr12-variable-sites.txt', 'r') as infile:
    sites = [int(x.strip()) - 1 for x in infile.readlines()]


for i in sites:
    fasta[i] = 'X'

fasta = ''.join(fasta)

forwardGuidePattern = re.compile("[ACTG]{20}[ACTG]GG")
reverseGuidePattern = re.compile("CC[ACTG][ACTG]{20}")
forwardGuides = list(re.finditer(forwardGuidePattern, fasta))
reverseGuides = list(re.finditer(reverseGuidePattern, fasta))

for i in forwardGuides:
    start = i.span()[0]
    start += 1
    stop = (start + 20)
    strand = "+"
    seq = i.group(0)
    pam = seq[-3:]
    seq = seq[:20]
    print("\t".join([str(x) for x in [seq, pam, start, stop, strand]]))


for i in reverseGuides:
    start = i.span()[1]
    stop += 1
    start = stop - 20
    strand = "-"
    seq = revcomp(i.group(0))
    pam = seq[-3:]
    seq = seq[:20]
    print("\t".join([str(x) for x in [seq, pam, start, stop, strand]]))