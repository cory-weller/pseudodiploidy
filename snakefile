import pandas as pd

samples = pd.read_table("samples.tsv").set_index("sample", drop=False)


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
        expand("data/input/{date}-{sample}-assigned/report.txt", date=DATE, sample=SAMPLES)
        #expand("data/input/{date}-{sample}.assembled.fastq", date=DATE, sample=SAMPLES)

rule get_reads:
    output: 
        "data/input/{date}-{sample}-R1-{stage}.fastq.gz",
        "data/input/{date}-{sample}-R2-{stage}.fastq.gz"
    threads: 1
    resources:
        mem_mb = 1024*1,
        runtime_min = 5
    params:
        filename1 = lambda wc: samples.loc[wc.sample, 'filename1'],
        filename2 = lambda wc: samples.loc[wc.sample, 'filename2'],
        sample = lambda wc: samples.loc[wc.sample, 'sample'],
        ida = lambda wc: samples.loc[wc.sample, 'ida'],
        idb1 = lambda wc: samples.loc[wc.sample, 'idb1'],
        idc1 = lambda wc: samples.loc[wc.sample, 'idc1'],
        idb2 = lambda wc: samples.loc[wc.sample, 'idb2'],
        idc2 = lambda wc: samples.loc[wc.sample, 'idc2']
    shell:
        """
        wget -O data/input/{params.filename1} https://onedrive.live.com/download?cid={params.ida}\\&resid={params.ida}%{params.idb1}\&authkey={params.idc1}
        wget -O data/input/{params.filename2} https://onedrive.live.com/download?cid={params.ida}\\&resid={params.ida}%{params.idb2}\&authkey={params.idc2}
        """

# snakemake -nr -j4 --latency-wait=180 --cluster="sbatch -c {threads} --mem={resources.mem_mb} --time={resources.runtime_min}"

rule pair_with_pear:
    input:
        r1="data/input/{date}-{sample}-R1-raw.fastq.gz",
        r2="data/input/{date}-{sample}-R2-raw.fastq.gz"
    output: "data/input/{date}-{sample}.assembled.fastq"
    threads: 2
    resources:
        mem_mb = 1024*4,
        runtime_min = 60*2
    container: "library://wellerca/pseudodiploidy/mapping:latest"
    shell:
        """
        pear -f {input.r1} -r {input.r2} -o data/input/{wildcards.date}-{wildcards.sample}
        rm data/input/{wildcards.date}-{wildcards.sample}.unassembled.{{forward,reverse}}.fastq
        rm data/input/{wildcards.date}-{wildcards.sample}.discarded.fastq
        """

rule assess_reads:
    input: "data/input/{date}-{sample}.assembled.fastq"
    #output: touch(".snakemake/assess_reads.done")
    output: "data/input/{date}-{sample}-assigned/report.txt"
    threads: 1
    resources:
        mem_mb = 1024*3,
        runtime_min = 10
    script: "src/test.py"