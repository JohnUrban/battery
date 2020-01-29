#!/bin/bash

mkdir $DIR
cd $DIR

PATH=~/software/miniasm/minitools/:~/software/miniasm/miniasm/:~/software/miniasm/minimap/:$PATH
PATH=~/software/miniasm/minitools/:~/software/miniasm/miniasm/:~/software/miniasm/minimap/:$PATH
PATH=~/software/racon/racon/bin/:$PATH

echo Overlapping; date
minimap -Sw5 -L100 -m0 -t$T $READS $READS | gzip -1 > reads.paf.gz

date; echo; echo; echo; echo
echo Assembly; date

miniasm -f $READS reads.paf.gz > reads.gfa
date; echo; echo; echo; echo

gfa2fasta reads.gfa > asm.fasta

minimap asm.fasta $READS > overlaps-for-racon.paf  
racon -t $T $READS overlaps-for-racon.paf asm.fasta asm.racon.fasta  

