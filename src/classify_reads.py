#!/usr/bin/env python

import sys
import time
import logging
import yaml
import regex as re
from yaml.loader import SafeLoader

logger = None

def yaml_load(yamlfilename):
    '''Loads a .yaml file and returns the contents as a dictionary'''
    with open(yamlfilename, 'r', encoding='utf-8') as yamlfile:
        contents = yaml.load(yamlfile, Loader=SafeLoader)
        return contents

def iterate_fasta(filename):
    '''Generator for reading through fasta, returning one entry at a time'''
    with open(filename, 'r', encoding='utf-8') as file:
        header = 'wt'
        seq = ''
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                yield (header, seq)
                header = line.strip(">")
                seq = ''
            else:
                seq += line
        yield (header, seq)

def classify_read():
    #thing

def main():
    '''Main'''
    startTime = time.time()
    config = yaml_load('config.yaml')
    # Initialize logging

    fasta_aln_infile = sys.argv[1]
    
    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level=config['logging'])
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(console_handler)
    run_time = (time.time() - startTime)
    logger.info("Total runtime: %s", run_time)

    for dna_read in iterate_fasta(fasta_aln_infile):
        print(dna_read)
    sys.exit()



if __name__ == '__main__':
    main()
