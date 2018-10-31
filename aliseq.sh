#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

snakemake --cluster-config $DIR/cluster.json \
          --jobscript $DIR/custombash.sh \
          --use-conda \
          --cluster "qsub -cwd -pe sharedmem {cluster.core} -l h_rt={cluster.time} -l h_vmem={cluster.vmem}" \
          --jobs 1000 \
          --snakefile $DIR/workflow.snk

