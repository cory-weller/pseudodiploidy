#!/usr/bin/env bash

snakemake -j 8 --latency-wait=180 --cluster="sbatch -c {threads} --mem={resources.mem_mb} --time={resources.runtime_min}"
