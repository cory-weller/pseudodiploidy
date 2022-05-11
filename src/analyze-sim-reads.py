#!/usr/bin/env python


import itertools
import regex as re
import gzip
import bisect
import sys
import os
import logging

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

def splitSeq(seq, geneName, regexLeft, regexRight):
    # Trim left constant sequence
    try:
        matchL = re.search(regexLeft, seq).captures()[0]
        beforeMatchL, afterMatchL = seq.split(matchL)
        matchR = re.search(regexRight, afterMatchL).captures()[0]
        beforeMatchR, afterMatchR = afterMatchL.split(matchR)
    except AttributeError:
        return(None)
    return([geneName, seq, beforeMatchL, matchL, beforeMatchR, matchR, afterMatchR])


import src.functions as f
CONFIG = f.yamlLoad('config.yaml')

LOGGING_LEVEL = CONFIG['logging']
DEBUG_ARGS = CONFIG['debug_args']
MATCH_LENGTH = CONFIG['match_length']
MISMATCH_TOLERANCE = CONFIG['mismatch_tolerance']
BUFFER_LENGTH = CONFIG['buffer_length']



args = sys.argv
if len(args) == 1:       # if no args supplied, e.g. debuging
    print('\n'.join(['running with debug arguments:'] + DEBUG_ARGS))
    args = DEBUG_ARGS



# set logging console handler
logger = logging.getLogger()
logger.setLevel(logging.NOTSET)
console_handler = logging.StreamHandler()
console_handler.setLevel(level=LOGGING_LEVEL)
console_handler_format = '%(asctime)s | %(filename)s %(levelname)s: %(message)s'
console_handler.setFormatter(logging.Formatter(console_handler_format))
logger.addHandler(console_handler)

# CEWID



# CRITICAL
# ERROR
# WARNING
# INFO
# DEBUG
#
# Import Arguments 
#
infile = args[1]        # infile is: /path/to/file/BATCH_SAMPLENAME.assembled.fastq


filestem = os.path.basename(infile).split('.assembled.fastq')[0] # filestem is: BATCH_SAMPLENAME.assembled.fastq
batch, samplename = filestem.split('-')

logger.debug('test')
logger.info('test')
logger.warning('test')
logger.error('test')
logger.critical('test')

genes = list(CONFIG['amplicons'].keys())

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





logger.info("Importing fastq sequences from %s" % infile)

goodSeqFile = "data/processed/%s-%s.good.txt" % (batch, samplename)
logger.info("Writing seqs with a match to %s" % goodSeqFile)

badSeqFile = "data/processed/%s-%s.bad.txt" % (batch, samplename)
logger.info("Writing seqs with no matches to %s" % badSeqFile)


with open(infile, 'r') as file, open(goodSeqFile, 'w') as goodOut, open(badSeqFile, 'w') as badOut:
    for lines in itertools.zip_longest(*[file]*4):
            break_out_flag = False
            rawSeq = lines[1].strip()   
            for gene in genes:
                result = splitSeq(rawSeq, gene, CONFIG['amplicons'][gene]['upRegexPattern'], CONFIG['amplicons'][gene]['downRegexPattern'])
                if result != None:
                    goodOut.write('\t'.join(result) + '\n')
                    break_out_flag = True
            if break_out_flag == True:
                continue
            badOut.write(rawSeq + '\n')


sys.exit()
