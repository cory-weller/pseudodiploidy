import pandas as pd
samples = pd.read_table("samples.tsv").set_index("sample", drop=False)
configfile: "config.yaml"

project_path = config['project_path']


SAMPLES = list(samples['sample'])   # SAMPLES: list of samples of given BATCH that will be processed
BATCH = config['run']['batch']      # BATCH: identifies the group of samples we're concerned with, e.g. date of sequencing run

# rule targets identifies the list of final output files to generate
rule targets:
    input:
        expand("data/processed/{batch}-{sample}.{types}.txt", batch=BATCH, sample=SAMPLES, types=['good','bad'])
        #expand("data/input/{batch}-{sample}.assembled.fastq", batch=BATCH, sample=SAMPLES)

rule get_reads:
    output: 
        "data/input/{batch}-{sample}-R1-{stage}.fastq.gz",
        "data/input/{batch}-{sample}-R2-{stage}.fastq.gz"
    threads: 1
    resources:
        mem_mb = 1024*1,
        runtime_min = 5
    params:
        filename1 = lambda wc: samples.loc[wc.sample, 'filename1'],
        filename2 = lambda wc: samples.loc[wc.sample, 'filename2'],
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
        r1="data/input/{batch}-{sample}-R1-raw.fastq.gz",
        r2="data/input/{batch}-{sample}-R2-raw.fastq.gz"
    output: "data/input/{batch}-{sample}.assembled.fastq"
    threads: 2
    resources:
        mem_mb = 1024*4,
        runtime_min = 60*2
    container: "library://wellerca/pseudodiploidy/mapping:latest"
    shell:
        """
        pear -f {input.r1} -r {input.r2} -o data/input/{wildcards.batch}-{wildcards.sample}
        rm data/input/{wildcards.batch}-{wildcards.sample}.unassembled.{{forward,reverse}}.fastq
        rm data/input/{wildcards.batch}-{wildcards.sample}.discarded.fastq
        """

rule assess_reads:
    input: "data/input/{batch}-{sample}.assembled.fastq"
    output: 
        o1="data/processed/{batch}-{sample}.good.txt",
        o2="data/processed/{batch}-{sample}.bad.txt"
    threads: 1
    resources:
        mem_mb = 1024*3,
        runtime_min = 10
    shell:
        """
        python3 src/analyze-sim-reads.py {input}
        """