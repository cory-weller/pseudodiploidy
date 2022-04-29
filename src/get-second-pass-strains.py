#!/usr/bin/env python

with open('data/input/strains-to-start-with.txt', 'r') as infile:
    startStrains = [x.strip() for x in infile.readlines()]

with open('data/input/first-pass-strains.txt', 'r') as infile:
    usedStrains = [x.strip() for x in infile.readlines()]

strainsForSecondPass = [x for x in startStrains if x not in usedStrains]

with open('data/input/second-pass-strains.txt', 'w') as outfile:
    outfile.write('\n'.join(strainsForSecondPass))