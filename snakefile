configfile: "config.yaml"
BOX = config['boxname'] + ':' + config['boxdir']
SAMPLES = config['samples']
PAIRS = config['pairs']
THREADS = config['threads']

rule all:
    input:
        #expand("data/input/{sample}_{pair}.fastq.gz", sample=SAMPLES, pair=PAIRS)
        expand("data/input/{sample}.assembled.fastq", sample=SAMPLES, pair=PAIRS)

rule get_reads:
    output:
        "data/input/{sample}_{pair}.fastq.gz"
    shell:
        "rclone copyto --progress {BOX}/{wildcards.sample}_L001_{wildcards.pair}_001.fastq.gz {output}"


rule pair_with_pear:
    input:
        r1="data/input/{sample}_R1.fastq.gz",
        r2="data/input/{sample}_R2.fastq.gz"
    output:
        "data/input/{sample}.assembled.fastq"
    shell:
        """
        module load pear
        pear -f {input.r1} -r {input.r2} -o data/input/{wildcards.sample}
        """
