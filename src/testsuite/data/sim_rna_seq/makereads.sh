#!/bin/bash

## everything will be uniformly expressed here....

G=../ecolitranscriptome/Escherichia_coli_k_12.ASM80076v1.cdna.all.fa
L=100
C=130
EM=0.001
EI=0.0001
ED=0.0001
MU=280
SD=60


O=rna1
fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -pe $MU $SD
O=rna2
fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -pe $MU $SD
rm *.i.fastq *.c.fastq

awk '{sub(/#1/,""); print}' rna1.1.fastq  > rna1.1
awk '{sub(/#2/,""); print}' rna1.2.fastq  > rna1.2
mv rna1.1 rna1.1.fastq
mv rna1.2 rna1.2.fastq

awk '{sub(/#1/,""); print}' rna2.1.fastq  > rna2.1
awk '{sub(/#2/,""); print}' rna2.2.fastq  > rna2.2
mv rna2.1 rna2.1.fastq
mv rna2.2 rna2.2.fastq


for i in 1 2; do r1=${PWD}/rna${i}.1.fastq; r2=${PWD}/rna${i}.2.fastq; echo -e ${r1}"\t"${r2}; done > reads.fofn
