#!/bin/bash

IN=../02-snpdecoligenome/02-ecoli-snpd-genome.fasta
OUT=03-ecoli.mutated

if [ ! -f $IN.fai ]; then samtools faidx $IN; fi

../../../third_party/SVsim/SVsim -i sv.commands -r $IN -o $OUT -W
