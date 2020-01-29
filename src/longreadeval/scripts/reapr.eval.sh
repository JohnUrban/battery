#!/bin/bash

## Need: REF, BASE, AGGRESSIVE (t/f), BAM


echo facheck
reapr facheck $REF ${BASE}_renamed

if $AGGRESSIVE ; then A="-break a=1"; else A=""; fi

echo pipeline

reapr pipeline ${BASE}_renamed.fa $BAM output_directory -stats s=1 -score l=$WINLEN -score f=$MININNER -score u=$FCDWINLEN -fcdrate l=$FCDWINLEN -fcdrate w=$MAXFCDSAMPLE ${A}


## do stuff
cd output_directory/
gunzip -c 03.score.per_base.gz | cut -f 3 | awkMean > per-base-mean-score.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py > broken_assembly.sizestats.txt
faSize -detailed 04.break.broken_assembly.fa | asm-stats.py -t > broken_assembly.sizestats.csv
cd ../

if $CLEAN; then bash ${SCRIPTS}/reapr.clean.sh ; fi
