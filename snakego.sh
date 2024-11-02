#!/bin/bash
#snakego.sh
snakefile=methylation_BS_EMseq.snake
key=$2

mkdir -p cluster cluster_log logs

echo "start at `date`"

snakemake --nolock \
	-s ${snakefile} -j 100 --latency-wait 1000 --force-use-threads  -pr --rerun-triggers mtime --rerun-incomplete \
	> cluster_log/0.${snakefile}.log 2> cluster/0.${snakefile}.err

echo "end at `date`"
