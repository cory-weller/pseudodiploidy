#!/usr/bin/env python

import sys
import gzip


vcf_filename = sys.argv[1]

previous_range = (0,1)


def overlapping(range1, range2):
    min_2 = min(range2)
    if min_2 <= max(range1) and min_2 >= min(range1):
        return True

start = 171320
end = 172320
OVERLAP = False
all_overlaps = []
current_overlaps = []
with gzip.open(vcf_filename, 'r') as infile:
    for line in infile:
        chromosome, pos, _, ref = line.decode().strip().split()[:4]
        pos = int(pos)
        if pos < start or pos > end:
            continue
        current_range = (pos, pos + len(ref) - 1)
        #print(pos, ref, current_range)
        #print(previous_range, current_range)
        if overlapping(previous_range, current_range):
            #OVERLAP = True
            #previous_range = (min(previous_range), max(current_range))
            current_overlaps.append(pos)
        else:           # Current range does not overlap
            # If current values being held are overlaps
            if len(current_overlaps) > 1:
                # append values
                all_overlaps.append(current_overlaps)
            previous_range = current_range
            current_overlaps = [pos]

for i in all_overlaps:
    print(i)

   


class variant:
    def __init__(self, vcf_row_text):
        self.chr, self.pos, _, self.ref, self.alts = vcf_row_text.strip().split()[:5]
        self.pos = int(self.pos)
        self.alts = self.alts.split(',')
        self.genos = [self.ref] + self.alts
        self.end = self.pos + len(self.ref) - 1

variant1 = variant("chromosome12    171934  .       CA      C,TA   asdfadsfasdf")
variant2 = variant("chromosome12    171935  .       A       T   asdfasdfasdfafsf")

def overlapping(v1, v2):
    if v1.chr == v2.chr & v1.end >= v2.pos:
        return True

def merge_rows(v1, v2):



var1 = variant('chr12',