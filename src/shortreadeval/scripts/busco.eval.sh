#!/bin/bash

## for use with busco.loop.slurm.in or eval script
## Need: TARGET, OUT, CPU

echo 1 $TARGET
echo 2 $OUT
echo 3 $CPU
echo 4 $BV1_LINEAGE
echo 5 $BV1_MODE

BUSCO_v1.22.py -in $TARGET -o $OUT -l $BV1_LINEAGE -m $BV1_MODE --cpu $CPU 



