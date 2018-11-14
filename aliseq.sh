#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

snakemake --cluster-config $DIR/scripts/cluster.json \
          --jobscript $DIR/scripts/custombash.sh \
          --use-conda \
          --cluster "qsub -cwd -pe sharedmem {cluster.core} -l h_rt={cluster.time} -l h_vmem={cluster.vmem}" \
          --jobs 1000 \
          --snakefile $DIR/scripts/workflow.snk

snakemake --report aliseq.report.html \
          --snakefile $DIR/scripts/workflow.snk 

# create the following directory if it doesn't exist already (-p prevents overwriting)
mkdir -p logs/jobreports

mv snakejob* logs/jobreports

