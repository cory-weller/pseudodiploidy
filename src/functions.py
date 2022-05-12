


def cleanSeq(seq):
    seq = seq.upper()
    for x in [' ', '\n', '\t']:
        seq = seq.replace(x, '')
    return(seq)


def reverseComplement(seq):
    pairs = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C",
        "N" : "N"
    }
    seq = seq.upper()
    for x in [' ', '\n', '\t']:
        seq = seq.replace(x, '')
    return(seq)
