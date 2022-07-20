#!/usr/bin/env python


# sed -r 's/>.*genomic\] \[[a-z]+=([^]]+?)].*$/>S288C_chr\1/g' ./data/external/S288C_reference_sequence_R64-3-1_20210421.fsa | src/testguides.py

import fileinput
import regex as re

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
    def __init__(self, match, strand, chromosome):
        self.strand = strand
        self.chromosome = chromosome
        if strand == '+':
            self.start = match.span()[0] + 1
            self.stop = (self.start + 19)
            self.seq = match.group(0)
            self.pam = self.seq[-3:]
            self.seq = self.seq[:20]
        elif strand == '-':
            self.start = match.span()[1]
            self.stop = self.start - 19
            self.seq = revcomp(match.group(0))
            self.pam = self.seq[-3:]
            self.seq = self.seq[:20]

def findGuides(chromosome, seq, forwardGuidePattern, reverseGuidePattern):
    for forwardGuide in re.finditer(forwardGuidePattern, seq):
        yield guide(forwardGuide, '+', chromosome)
    for reverseGuide in re.finditer(reverseGuidePattern, seq):
        yield guide(reverseGuide, '-', chromosome)

forwardGuidePattern = re.compile("[ACTG]{20}[ACTG]GG")
reverseGuidePattern = re.compile("CC[ACTG][ACTG]{20}")

text = ''.join(fileinput.input()).split('>')[1:]


for entry in text:
    header, seq = entry.split('\n', 1)
    seq = seq.replace('\n', '')
    guides = findGuides(header, seq, forwardGuidePattern, reverseGuidePattern)
    for i in guides:
        print(i.seq, i.strand, i.pam, i.chromosome, i.start, i.stop, sep='\t')

