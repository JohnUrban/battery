#!/bin/bash
###########################

mkdir $DIR
cd $DIR

module load blasr/2015Oct22-8cc8621 
export PATH=~/software/abruijn/ABruijn/:~/software/abruijn/ABruijn/bin:$PATH

OUTDIR=asm

abruijn.py $READS $OUTDIR $COV -t ${T} --iterations 1 ##--resume


