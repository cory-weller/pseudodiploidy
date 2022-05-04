#!/usr/bin/env python

import itertools
import regex as re
import gzip
import bisect
import sys

#args = sys.argv
args = ['analyze-reads.py', 'barcodes.txt']
infile = args[1]
#outfile = sys.argv[2]

class fastq_read:
    def __init__(self, lines):
        self.id, self.seq, self.line3, self.quality = lines

class read:
    def __init__(self, name, rawSeq):
        self.name = name
        self.rawSeq = rawSeq
        self.upBarcode = ''
        self.downBarcode = ''
        self.repairTemplate = ''
        self.count = 1
    #@property
    #def increment(self):
    #    self.count += 1
    def extractSeq(self, fullPattern, upstreamPattern, downstreamPattern):
        self.matches = re.search(fullPattern, self.rawSeq)
        if self.matches == None:
            self.RTseq = None
            #print("no matches")
        elif len(self.matches.captures()) > 1:
            self.RTseq = None
            #print("more than one pattern match, do what now?")
        else:
            self.extractedSeq = self.matches.captures()[0]
            self.preMatches = re.search(upstreamPattern, self.extractedSeq)
            self.preMatchString = self.preMatches.captures()[0]
            self.preMatchEnd = self.preMatches.ends()[0]
            self.postMatches = re.search(downstreamPattern, self.extractedSeq)
            self.postMatchString = self.postMatches.captures()[-1]
            self.postMatchStart = self.postMatches.starts()[-1]
            self.RTseq = self.extractedSeq[self.preMatchEnd : self.postMatchStart]
    def printRT(self):
        print(self.RTseq)


class reads:
    def __init__(self, f):
        self.all = []
        if f.endswith('.gz'):
            self.gzip = True
            self.filename = f.rstrip(".gz")
        else:
            self.gzip = False
            self.filename = f
        self.extension = self.filename.split('.')[-1].lower()
        if self.extension in ["fastq", "fq"]:
            self.readType = "fastq"
        elif self.extension in ["fasta", "fa"]:
            self.readType = "fasta"
        self.entries = {}
    def importGzipFastq(self):
        print("Importing gzipped fastq sequences from %s" % self.filename)
        with gzip.open(self.filename + '.gz', 'rb') as f:
            for lines in itertools.zip_longest(*[f]*4):      
                rawSeq = lines[1].strip().decode()
                try:
                    print("\t".join(splitSeq(i.rawSeq)))
                except:
                    pass
                #self.all.append(read('nullName', rawSeq))
    def importFastq(self):
        print("Importing fastq sequences from %s" % self.filename)
        with open(self.filename, 'r') as f:
            for lines in itertools.zip_longest(*[f]*4):      
                rawSeq = lines[1].strip()
                readName = lines[0].strip()
                self.all.append(read(readName, rawSeq))
    def importGzipFasta(self):
        print("Importing gzipped fasta sequences from %s" % self.filename)
        with gzip.open(self.filename, 'rb') as f:
            samples = f.read().split(">")[1:]
            samples = [x.strip().split("\n") for x in samples]
            samples = [(x[0], ''.join(x[1:])) for x in samples]
    def importFasta(self):
        print("Importing fasta sequences from %s" % self.filename)
        with open(self.filename, 'r') as f:
            samples = f.read().split(">")[1:]
            samples = [x.strip().split("\n") for x in samples]
            samples = [(x[0], ''.join(x[1:])) for x in samples]
            for i in samples:
                readName = i[0]
                rawSeq = i[1]
                self.all.append(read(readName, rawSeq))

plasmidLtolerance = 1
plasmidRtolerance = 1
sublibraryLtolerance = 1
sublibraryRtolerance = 2
sublibraryLmin = 15
sublibraryLmax = 15
sublibraryRmin = 15
sublibraryRmax = 15

plasmidL = ['ACTACCGGCTGATATCATC']
sublibraryL = ['CGGAGGCAGGAGGGT', 'GCGGGTCACTTGGGT', 'CCGGGACCAGGCTCT']
sublibraryR = ['GGACTAGCGGCCGTC', 'GACAGTGGACCGGGC', 'GTCCCGCCGATCCTG']
plasmidR = ['CCTGAGTAACCGGTTC']

pPlasmidL = re.compile(r"""\L<plasmidL>{d<=%s}""" % (plasmidLtolerance), plasmidL=plasmidL)
pPlasmidR = re.compile(r"""\L<plasmidR>{d<=%s}""" % (plasmidRtolerance), plasmidR=plasmidR)
pSublibraryL = re.compile(r"""\L<sublibraryL>""", sublibraryL=sublibraryL)
pSublibraryR = re.compile(r"""\L<sublibraryR>""", sublibraryR=sublibraryR)

def splitSeq(myseq):
    # Trim left plasmid sequence
    matchPlasmidL = re.search(pPlasmidL, myseq).captures()[0]
    myseq = ''.join(myseq.split(matchPlasmidL)[1:])
    # Trim right plasmid sequence
    matchPlasmidR = re.search(pPlasmidR, myseq).captures()[-1]
    myseq = ''.join(myseq.split(matchPlasmidR)[:-1])
    # Split on (sublibraryL) into (L barcode, remainder) 
    matchSublibraryL = re.search(pSublibraryL, myseq).captures()[0]
    barcodeL, myseq = myseq.split(matchSublibraryL)
    # Split on (sublibraryR) into (RT, barcodeR)
    matchSublibraryR = re.search(pSublibraryR, myseq).captures()[-1]
    RT, barcodeR = myseq.split(matchSublibraryR)
    return([matchPlasmidL, barcodeL, matchSublibraryL, RT, matchSublibraryR, barcodeR, matchPlasmidR])


print("Importing gzipped fastq sequences from %s" % infile)
with open(infile, 'r') as f:
    with open(outfile, 'w') as o:
        for lines in itertools.zip_longest(*[f]*4):      
            rawSeq = lines[1].strip()
            try:
                o.write("\t".join(splitSeq(rawSeq))+"\n")
            except:
                pass
                #print("""error with %s""" % (rawSeq))
