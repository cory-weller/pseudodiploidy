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


def match_sequence_fuzzy(seq, amplicons):
    '''returns array for whether or not the given seqs match any patterns'''
    for amplicon in amplicons:
        for pattern in amplicon.patterns:
            if re.search(pattern, seq):
                return amplicon.name
    return None

def assign_sequence(seq, match, amplicons):
    '''assigns sequencing read (seq) corresponding to amplicon (match) as wt or edit'''
    for amplicon in amplicons:
        if amplicon.name == match:
            return amplicon.classify_read(seq)

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
                max_errors,
                cut_site,
                dist_around_edit
    ):
        self.name = name
        self.upstream = upstream
        self.downstream = downstream
        self.wt = wt
        self.edit = edit
        self.dist_around_edit = dist_around_edit
        self.seq = self.upstream[pad_left : -pad_right]
        self.patterns = [re.compile(fr"({self.seq}){{e<={error}}}") for error in max_errors]
        self.cut_site = cut_site
        self.slice = slice(self.cut_site - self.dist_around_edit, 
                           self.cut_site + self.dist_around_edit)
        self.wt_seq = (self.upstream + self.wt + self.downstream)[self.slice]
        self.edit_seq = (self.upstream + self.edit + self.downstream)[self.slice]
    def classify_read(self, read):
        if self.edit_seq in read:
            return 'edit'
        if self.wt_seq in read:
            return 'wt'
        else:
            return 'other'


def main():
    '''Main'''
    start_time = time.time()
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
        amplicon = Amplicon(name = gene,
                    upstream = config['amplicons'][gene]['upstream'],
                    downstream = config['amplicons'][gene]['downstream'],
                    wt = config['amplicons'][gene]['wt'],
                    edit = config['amplicons'][gene]['edit'],
                    pad_left = config['pad_left'],
                    pad_right = config['pad_right'],
                    max_errors={int(len(config['amplicons'][gene]['upstream'])*i) for i in mismatch_tolerances},
                    cut_site = config['amplicons'][gene]['cut_site'],
                    dist_around_edit = config['dist_around_edit']
        )
        amplicons.append(amplicon)

    for amplicon in amplicons:
        logger.info("adding amplicon: %s", amplicon.name)
        logger.info("adding perfect match pattern: %s", amplicon.patterns)
        logger.info("distance around cut site considered: %s", amplicon.dist_around_edit)
        logger.info("sequence diagnostic of wild type: %s", amplicon.wt_seq)
        logger.info("sequence diagnostic of edit: %s", amplicon.edit_seq)

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
            if match is None:
                match = 'nomatch'
                read_classification = 'other'
            else:
                read_classification = assign_sequence(seq, match, amplicons)
            print(f"{count}\t{read_classification}\t{seq}", file=files[match])
    run_time = (time.time() - start_time)
    logger.info("Total runtime: %s", run_time)
    sys.exit()


if __name__ == '__main__':
    main()
