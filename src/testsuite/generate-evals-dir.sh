#!/bin/bash
set -e

	
##################################
### EVALS ########################
##################################

mkdir -p evals
cp util/setup.cfg evals/
cp data/sim_rna_seq/reads.fofn evals/rnaseqreads.fofn

while read relpath; do
 readlink -f $relpath 
done < util/relative-input.fofn > evals/input.fofn

while read relpath; do
 readlink -f $relpath 
done < util/relative-input.extended.fofn > evals/input.extended.fofn

