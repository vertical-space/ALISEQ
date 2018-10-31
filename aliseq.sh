#!/bin/bash

snakemake --cluster-config cluster.json \
          --jobscript custombash.sh \
          --use-conda \
          --cluster "qsub -cwd -pe sharedmem {cluster.core} -l h_rt={cluster.time} -l h_vmem={cluster.vmem}" \
          --jobs 1000 \
          --snakefile workflow.snk

