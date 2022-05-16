#!/usr/bin/env python
'''docstring'''
import sys
import os
import logging
import contextlib
import itertools
import yaml
import regex as re
from yaml.loader import SafeLoader

#import inspect
#from Bio import AlignIO

def yaml_load(yamlfilename):
    '''Loads a .yaml file and returns the contents as a dictionary'''
    with open(yamlfilename, 'r', encoding='utf-8') as yamlfile:
        contents = yaml.load(yamlfile, Loader=SafeLoader)
        return contents


def clean_seq(seq):
    '''removes whitespace from sequences and converts to upper-case'''
    seq = seq.upper()
    for i in [' ', '\n', '\t']:
        seq = seq.replace(i, '')
    return seq


def reverse_complement(seq):
    '''returns reverse complement of a sequence (after converting to uppercase
    and removing whitespace. N nucleotides stay as N.'''
    pairs = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C",
        "N" : "N"
    }
    seq = seq.upper()
    for i in [' ', '\n', '\t']:
        seq = seq.replace(i, '')
    try:
        reverse_seq = ''.join([pairs[i] for i in seq[::-1]])
    except KeyError as err:
        logger.error('The sequence being reverse complemented,'
                     ' %s contains characters other than [ACTGN]', seq)
        sys.exit(f"ERROR: {err}")
    return reverse_seq


def split_seq(seq, geneName, regexLeft, regexRight):
    ''' Trims sequence before and after defined constant region regex patterns'''
    # Trim left constant sequence
    try:
        match_L = re.search(regexLeft, seq).captures()[0]
        before_match_L, after_match_L = seq.split(match_L)
        match_R = re.search(regexRight, after_match_L).captures()[0]
        before_match_R, after_match_R = after_match_L.split(match_R)
    except AttributeError:
        return None
    return([geneName, seq, before_match_L, match_L, before_match_R, match_R, after_match_R])


def check_pythonpath():
    '''Confirms current directory is in $PYTHONPATH'''
    try:
        assert os.environ['PYTHONPATH'] == os.getcwd()
    except AssertionError:
        sys.exit(
        f'ERROR while running {os.path.basename(sys.argv[0])}:'
            ' $PYTHONPATH environmental variable does not match current directory.'
            ' Ensure you are in the top level of the project directory,'
            ' and that $PYTHONPATH is set to the top level of the project directory.'
            ' e.g., export PYTHONPATH=/path/to/pseudodiploidy/'
        )


def match_sequence_perfect(seq, amplicons):
    '''returns array for whether or not the given seqs match any patterns'''
    output = []
    for i in amplicons:
        if re.search(i.regex_pattern1, seq):
            output.append(i.name)
    if len(output) > 0:
        return output
    return None

def match_sequence_fuzzy(seq, amplicons):
    '''returns array for whether or not the given seqs match any patterns'''
    output = []
    for i in amplicons:
        if re.search(i.regex_pattern2, seq):
            return(i.name)
    if len(output) > 0:
        return output
    return ['nomatch']


class Amplicon:
    def __init__(self,
                name,
                upstream,
                downstream,
                wt,
                edit,
                pad_left,
                pad_right,
                max_errors
    ):
        self.name = name
        self.upstream = upstream
        self.downstream = downstream
        self.wt = wt
        self.edit = edit
        self.seq = self.upstream[pad_left : -pad_right]
        self.regex_pattern1 = re.compile(fr"({self.seq})")
        self.regex_pattern2 = re.compile(fr"({self.seq}){{e<={max_errors}}}")

if __name__ == '__main__':
    # Load config and 

    CONFIG = yaml_load('config.yaml')
    DEBUG_ARGS = CONFIG['debug_args']

    # Initialize logging
    LOGGING_LEVEL = CONFIG['logging']
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level=LOGGING_LEVEL)
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(console_handler)

    # Parse args
    myargs = sys.argv
    if len(myargs) == 1:       # if no args supplied, e.g. debuging
        print('\n'.join(['running with debug arguments:'] + DEBUG_ARGS))
        myargs = DEBUG_ARGS

    fastq_name = myargs[1]        # /path/to/file/BATCH_SAMPLENAME.assembled.fastq
    logger.info("fastq input filename: %s", fastq_name)
    # filestem formatted/parsed with hyphen delimiter, as BATCH-SAMPLE.assembled.fastq
    filestem = os.path.basename(fastq_name).split('.assembled.fastq')[0] 
    logger.info("filestem: %s", filestem)
    batch, samplename = filestem.split('-')
    logger.info("batch: %s", batch)
    logger.info("samplename: %s", samplename)

    out_folder = '/'.join(fastq_name.split('/')[:-1]) + '/'

    # Initialize amplicons
    MISMATCH_TOLERANCE = float(CONFIG['pct_mismatch_tolerance'])
    gene_names = list(CONFIG['amplicons'].keys())
    amplicons = []
    for gene in gene_names:
        amplicon = Amplicon(gene,
                    CONFIG['amplicons'][gene]['upstream'],
                    CONFIG['amplicons'][gene]['downstream'],
                    CONFIG['amplicons'][gene]['wt'],
                    CONFIG['amplicons'][gene]['edit'],
                    CONFIG['pad_left'],
                    CONFIG['pad_right'],
                    int(len(CONFIG['amplicons'][gene]['upstream'])* MISMATCH_TOLERANCE)
        )
        amplicons.append(amplicon)
    
    for amplicon in amplicons:
        logger.info("adding amplicon: %s", amplicon.name)
        logger.info("adding perfect match pattern: %s", amplicon.regex_pattern1)
        logger.info("adding fuzzy match pattern: %s", amplicon.regex_pattern2)
    
    # Evaluate reads

    with contextlib.ExitStack() as stack, \
    open(fastq_name, 'r', encoding='utf-8') as infile:
        files = {}
        unmatched = []
        for gene in gene_names +  ['nomatch']:
            filename = f'{samplename}-{gene}.txt'
            files[gene] = stack.enter_context(open(filename, 'w', encoding='utf-8'))
        # Iterate over fastq lines for perfect matches
        logger.info("iterating over raw reads")
        for line in itertools.zip_longest(*[infile]*4):
            seq = line[1].strip()
            matches = match_sequence_perfect(seq, amplicons)
            if matches is None:
                unmatched.append(seq)
            else:
                for match in matches:
                    print(seq, file=files[match])
        # Iterate over unmatched sequences with fuzzy tolerance
        logger.info("iterating over %s unmatched reads with fuzzy tolerance", len(unmatched))
        for seq in unmatched:
            matches = match_sequence_fuzzy(seq, amplicons)
            for match in matches:
                print(match, file=files[match])

    sys.exit()