#!/usr/bin/env python
'''docstring'''
import sys
import os
import logging
import contextlib
import itertools
# import gzip
# import bisect
import yaml
import src.functions as f
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
    except KeyError as error:
        logger.error('sequence %s being reverse complemented contains something outside of [ACTGN]', seq)
        sys.exit(f"ERROR: {error}")
    return reverse_seq


CONFIG = yaml_load('config.yaml')
LOGGING_LEVEL = CONFIG['logging']

# set logging console handler
logger = logging.getLogger()
logger.setLevel(logging.NOTSET)
console_handler = logging.StreamHandler()
console_handler.setLevel(level=LOGGING_LEVEL)
console_handler_format = '%(levelname)s: %(message)s'
console_handler.setFormatter(logging.Formatter(console_handler_format))
logger.addHandler(console_handler)

# set logging shutdown handler

# try:
#     logger.info("Attempting to load src/functions.py")
#     import src.functions as f
# except ModuleNotFoundError:
#     logger.warning("src/functions.py not found in $PYTHONPATH.")
#     logger.warning(f"Setting os.environ['PYTHONPATH'] to current directory {os.getcwd()} and trying again." % (os.getcwd()))
#     try:
#         os.environ['PYTHONPATH'] = os.getcwd()
#         import src.functions as f
#     except ModuleNotFoundError as e:
#         logger.error( "Could not load src/functions.py. Confirm you are running from the top level of this project directory and/or set $PYTHONPATH to the top level of this project directory and try again.")
#         sys.exit("ModuleNotFoundError: " + str(e))
#     else:
#         logger.info("Successfully loaded src/functions.py after automatically setting $PYTHONPATH")
# else:
#     logger.info("Successfully loaded src/functions.py")


try:
    assert os.environ['PYTHONPATH'] == os.getcwd()
except AssertionError as err:
    sys.exit("""
    ERROR while running %s:
    $PYTHONPATH environmental variable does not match current directory. Ensure you are in the top level of the project directory,
    and that $PYTHONPATH is set to the top level of the project directory.
    e.g., export PYTHONPATH=/path/to/pseudodiploidy/

    $PYTHONPATH: %s
    current directory: %s

    """ % (os.path.basename(sys.argv[0]), os.environ['PYTHONPATH'], os.getcwd(), ))

def split_seq(seq, geneName, regexLeft, regexRight):
    ''' Trims sequence before and after defined constant region regex patterns'''
    # Trim left constant sequence
    try:
        matchL = re.search(regexLeft, seq).captures()[0]
        beforeMatchL, afterMatchL = seq.split(matchL)
        matchR = re.search(regexRight, afterMatchL).captures()[0]
        beforeMatchR, afterMatchR = afterMatchL.split(matchR)
    except AttributeError:
        return None
    return([geneName, seq, beforeMatchL, matchL, beforeMatchR, matchR, afterMatchR])





MATCH_LENGTH = CONFIG['match_length']
MISMATCH_TOLERANCE = CONFIG['mismatch_tolerance']
BUFFER_LENGTH = CONFIG['buffer_length']
DEBUG_ARGS = CONFIG['debug_args']

args = sys.argv
if len(args) == 1:       # if no args supplied, e.g. debuging
    print('\n'.join(['running with debug arguments:'] + DEBUG_ARGS))
    args = DEBUG_ARGS


infile = args[1]        # infile is: /path/to/file/BATCH_SAMPLENAME.assembled.fastq


filestem = os.path.basename(infile).split('.assembled.fastq')[0] # filestem is: BATCH_SAMPLENAME.assembled.fastq
batch, samplename = filestem.split('-')

# logger.debug('test')
# logger.info('test')
# logger.warning('test')
# logger.error('test')
# logger.critical('test')

genes = list(CONFIG['amplicons'].keys())
genes.append('nomatch')

for gene in genes:
    logger.debug('Generating regex pattern for ' + gene)

    upNucleotides = CONFIG['amplicons'][gene]['upstream'][::-1][BUFFER_LENGTH:(BUFFER_LENGTH + MATCH_LENGTH)][::-1]       # reverse, take first MATCH_LENGTH chars, reverse again
    logger.debug('upNucleotides: ' + upNucleotides)

    upString = r"%s{s<=%s}" % (upNucleotides,MISMATCH_TOLERANCE)
    logger.debug('upString: ' + upString)

    upRegexPattern = re.compile(upString)
    CONFIG['amplicons'][gene]['upRegexPattern'] = upRegexPattern
    logger.debug('upRegexPattern: ' + str(upRegexPattern))

    downNucleotides = CONFIG['amplicons'][gene]['downstream'][BUFFER_LENGTH:(BUFFER_LENGTH + MATCH_LENGTH)]
    logger.debug('downNucleotides: ' + downNucleotides)

    downString = r"%s{s<=%s}" % (downNucleotides,MISMATCH_TOLERANCE)
    logger.debug('downString: ' + downString)

    downRegexPattern = re.compile(downString)
    CONFIG['amplicons'][gene]['downRegexPattern'] = downRegexPattern
    logger.debug('downRegexPattern: ' + str(downRegexPattern))

for gene in genes:
    logger.debug('Generating regex pattern for ' + gene)

    upNucleotides = CONFIG['amplicons'][gene]['upstream'][::-1][BUFFER_LENGTH:(BUFFER_LENGTH + MATCH_LENGTH)][::-1]       # reverse, take first MATCH_LENGTH chars, reverse again
    logger.debug('upNucleotides: ' + upNucleotides)

    downNucleotides = CONFIG['amplicons'][gene]['downstream'][BUFFER_LENGTH:(BUFFER_LENGTH + MATCH_LENGTH)]
    logger.debug(f'downNucleotides: {downNucleotides}')

    regexPattern = r"%s.*%s{s<=%s}" % (upNucleotides, downNucleotides, MISMATCH_TOLERANCE)

    CONFIG['amplicons'][gene]['downRegexPattern'] = downRegexPattern
    logger.debug(f'downRegexPattern: {downRegexPattern}')




logger.info("Importing fastq sequences from %s" , infile)

goodSeqFile = "data/processed/%s-%s.good.txt" % (batch, samplename)
logger.info("Writing seqs with a match to %s" , goodSeqFile)

badSeqFile = "data/processed/%s-%s.bad.txt" % (batch, samplename)
logger.info("Writing seqs with no matches to %s", badSeqFile)


with open(infile, 'r') as infile, open(goodSeqFile, 'w', encoding='utf-8') as goodOut, open(badSeqFile, 'w', encoding='utf-8') as badOut:
    for lines in itertools.zip_longest(*[infile]*4):
        breakoutFlag = False
        rawSeq = lines[1].strip() 
        for gene in genes:
            result = split_seq(rawSeq, gene, CONFIG['amplicons'][gene]['upRegexPattern'], CONFIG['amplicons'][gene]['downRegexPattern'])
            if result is not None:
                goodOut.write('\t'.join(result) + '\n')
                breakoutFlag = True
        if breakoutFlag:
            continue
        badOut.write(rawSeq + '\n')

genes = ['spo13','can1','his2']

genes.append('nomatch')

sample = 'S2'



with contextlib.ExitStack() as stack,  open('data/input/100-tiny.assembled.fastq', 'r') as infile:
    files = {}
    for gene in genes:
        filename = f'{sample}-{gene}.txt'
        files[gene] = stack.enter_context(open(filename, 'w'))
    # set up file handling
    for line in itertools.zip_longest(*[infile]*4):
        rawseq = line[1].strip()
        if 'TTGGATGTTAATTGTTCCGATTATGAG' in rawseq:
            match = 'spo13'
        elif 'CTTTACCTTGTCATAGCACACAC' in rawseq:
            match = 'his2'
        elif 'AGAAAACTGTGAAAGAGGATGTAACAGGGA' in rawseq:
            match = 'can1'
        else:
            match = 'nomatch'
        print(rawseq, file=files[match])

sys.exit()
