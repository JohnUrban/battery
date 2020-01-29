#!/bin/bash

## for use with busco.loop.slurm.in or eval script
## Need: FASTA, OUT, CPU

for VAR in FASTA OUT CPU LINEAGE BV3_MODE REGIONLIMIT PATH; do echo $VAR ${!VAR}; echo; done; echo

run_BUSCO.py --in $FASTA -o $OUT -l $LINEAGE -m $BV3_MODE --cpu $CPU --limit $REGIONLIMIT


