#!/bin/bash

# to run snakemake as batch job
module load snakemake || exit 1
module load bedtools/2.27.1 

mkdir -p 00log 
sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s /home/mcgaugheyd/git/ipsc_rpe_atac/Snakefile \
-pr --local-cores 8 --jobs 1999 \
--configfile $1 \
--cluster-config /home/mcgaugheyd/git/ipsc_rpe_atac/cluster.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 1 \
--resources parallel=4
