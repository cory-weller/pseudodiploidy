#!/usr/bin/env python

seqs = {}

import sys
clustal_filename = sys.argv[1]

with open(clustal_filename, 'r', encoding='utf-8') as clustal:
    for line in clustal:
        if "_" not in line:
            continue
        text = line.strip().split()
        strain = text[0].split('_')[0]
        seq = text[1]
        if strain not in seqs:
            seqs[strain] = []
        seqs[strain].append(seq)

for strain in list(seqs.keys()):
    concat_seq = ''.join(seqs[strain])
    print(f'{strain}\t{concat_seq}')
