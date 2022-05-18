#!/usr/bin/env python
'''docstring'''
import sys
import os
import logging
import contextlib
import itertools
import yaml
import time
import regex as re
from yaml.loader import SafeLoader

logger = None

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


# def match_sequence_perfect(seq, amplicons):
#     '''returns array for whether or not the given seqs match any patterns'''
#     output = []
#     for i in amplicons:
#         if re.search(i.regex_pattern1, seq):
#             output.append(i.name)
#     if len(output) > 0:
#         return output
#     return None

def match_sequence_fuzzy(seq, amplicons):
    '''returns array for whether or not the given seqs match any patterns'''
    for amplicon in amplicons:
        for pattern in amplicon.patterns:
            if re.search(pattern, seq):
                return amplicon.name
    return None


class Amplicon:
    '''describes information in config.yaml for identifying various sequenced amplicons'''
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
        self.patterns = [re.compile(fr"({self.seq}){{e<={error}}}") for error in max_errors]


def main():
    '''Main'''
    startTime = time.time()
    config = yaml_load('config.yaml')

    # Initialize logging
    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level=config['logging'])
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(console_handler)

    # Parse args
    myargs = sys.argv
    if len(myargs) == 1:       # if no args supplied, e.g. debuging
        myargs = config['debug_args']
        logger.debug("Running with debug arguments defined in config.yaml")
    logger.info("Arguments: %s", myargs)

    fastq_name = myargs[1]        # /path/to/file/BATCH_SAMPLENAME.assembled.fastq
    logger.info("fastq input filename: %s", fastq_name)
    # filestem formatted/parsed with hyphen delimiter, as BATCH-SAMPLE.assembled.fastq
    filestem = os.path.basename(fastq_name).split('.assembled.fastq')[0]
    logger.info("filestem: %s", filestem)
    batch, samplename = filestem.split('-')
    logger.info("batch: %s", batch)
    logger.info("samplename: %s", samplename)

    # Initialize amplicons
    mismatch_tolerances = config['pct_mismatch_tolerances']
    logger.info("Mismatch Tolerances: %s", mismatch_tolerances)
    gene_names = list(config['amplicons'].keys())
    amplicons = []
    for gene in gene_names:
        amplicon = Amplicon(gene,
                    config['amplicons'][gene]['upstream'],
                    config['amplicons'][gene]['downstream'],
                    config['amplicons'][gene]['wt'],
                    config['amplicons'][gene]['edit'],
                    config['pad_left'],
                    config['pad_right'],
                    {int(len(config['amplicons'][gene]['upstream'])*i) for i in mismatch_tolerances}
        )
        amplicons.append(amplicon)

    for amplicon in amplicons:
        logger.info("adding amplicon: %s", amplicon.name)
        logger.info("adding perfect match pattern: %s", amplicon.patterns)

    # Evaluate reads

    with contextlib.ExitStack() as stack, \
    open(fastq_name, 'r', encoding='utf-8') as infile:
        files = {}
        for gene in gene_names +  ['nomatch']:
            filename = f'data/processed/{batch}-{samplename}-{gene}.counts'
            files[gene] = stack.enter_context(open(filename, 'w', encoding='utf-8'))
        # Iterate over fastq lines for perfect matches
        logger.info("iterating over raw reads")
        i = 0
        # for line in itertools.zip_longest(*[infile]*4):  # for iterating over fastq
        for line in infile:     # for iterating over table of <count> <read>
            i += 1
            if(i % 10000) == 0:
                logger.info("%s reads done", i)
            count, seq = line.strip().split()
            match = match_sequence_fuzzy(seq, amplicons)
            # here
            if match is None:
                match = 'nomatch'
            print(f"{count}\t{seq}", file=files[match])
    run_time = (time.time() - startTime)
    logger.info("Total runtime: %s", run_time)
    sys.exit()


if __name__ == '__main__':
    main()
