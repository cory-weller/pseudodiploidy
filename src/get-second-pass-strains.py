#!/usr/bin/env python

import sys

starting_strain_file = sys.argv[1]
used_strains_file = sys.argv[2]

with open(starting_strain_file, 'r', encoding='utf-8') as infile:
    startStrains = [x.strip() for x in infile.readlines()]

with open(used_strains_file, 'r', encoding='utf-8') as infile:
    usedStrains = [x.strip() for x in infile.readlines()]

strainsForSecondPass = [x for x in startStrains if x not in usedStrains]

for i in strainsForSecondPass:
    print(i)