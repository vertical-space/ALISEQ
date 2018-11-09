#!/bin/bash

# create the following directories if they don't exist already (-p prevents overwriting)
mkdir -p stats logs/jobreports

# get the full path to the toplevel project directory (this is the toplevel directory pulled down from github, called ALISEQ)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

snakemake --cluster-config $DIR/scripts/cluster.json \
          --jobscript $DIR/scripts/custombash.sh \
          --use-conda \
          --cluster "qsub -cwd -pe sharedmem {cluster.core} -l h_rt={cluster.time} -l h_vmem={cluster.vmem}" \
          --jobs 1000 \
          --snakefile $DIR/scripts/simulate_reads.snk

snakemake --report simulate_reads.report.html \
          --snakefile $DIR/scripts/simulate_reads.snk 

mv snakejob* logs/jobreports

