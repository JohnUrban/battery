#!/bin/bash

G=../../asms/refs/01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta
O=ilmn
L=100
C=100
EM=0.001
EI=0.0001
ED=0.0001
MU=400
SD=60


## From Canu
fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -pe $MU $SD
rm *.i.fastq *.c.fastq
