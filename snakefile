configfile: "config.yaml"
BOX = config['boxname'] + ':' + config['boxdir']
SAMPLES = config['samples']
PAIRS = config['pairs']
THREADS = config['threads']

#import yaml
#from yaml.loader import SafeLoader
#with open('config.yaml', 'r') as infile:
#    config = yaml.load(infile, Loader=SafeLoader)


DATE = config['run']['date']
SAMPLES = [key for key in config['run']['sample'].keys()]
READS = [key for key in config['run']['reads']]
EXTENSION = config['run']['extension']
STAGE = config['run']['stage']
CID = config['run']['cid']

rule all:
    input:
        expand("data/input/{date}-{sample}.assembled.fastq", date=DATE, sample=SAMPLES)
    

rule retrieve_reads:
    output:
        R1="data/input/{date}-{sample}-R1-{stage}.fastq.gz",
        R2="data/input/{date}-{sample}-R2-{stage}.fastq.gz"
    params:
        fid1 = lambda wildcards: config['run']['sample'][wildcards.sample]['R1']['fid'],
        authkey1 = lambda wc: config['run']['sample'][wc.sample]['R1']['authkey'],
        fid2 = lambda wildcards: config['run']['sample'][wildcards.sample]['R2']['fid'],
        authkey2 = lambda wc: config['run']['sample'][wc.sample]['R2']['authkey']
    shell:
        "wget -O {output.R1} https://onedrive.live.com/download?cid={CID}\&resid={CID}%{params.fid1}\&authkey={params.authkey1}"
        "wget -O {output.R2} https://onedrive.live.com/download?cid={CID}\&resid={CID}%{params.fid2}\&authkey={params.authkey2}"

rule pair_with_pear:
    input:
        r1="data/input/{sample}-R1-raw.fastq.gz",
        r2="data/input/{sample}-R2-raw.fastq.gz"
    output:
        "data/input/{sample}.assembled.fastq"
    shell:
        """
        module load pear
        pear -f {input.r1} -r {input.r2} -o data/input/{wildcards.sample}
        """