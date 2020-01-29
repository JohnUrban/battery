#!/bin/bash

G=../../asms/refs/01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta
O=pacbio
L=8000
C=50
EM=0.05
EI=0.05
ED=0.05
#fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

LENS="18390 16199 14835 13825 13023 12324 11682 11035 10363 9681 8992 8312 7638 6958 6269 5551 4756 3808 2613 1000"
C=2.5
for L in $LENS; do 
  O=tmp-${L}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se
done 2> errfile
cat tmp*fastq > pacbio.s.fastq
rm tmp*fastq
fastq2faqual.py -f pacbio.s.fastq --fa --out pacbio.s

#$N5[1] 18390
#$N10[1] 16199
#$N15[1] 14835
#$N20[1] 13825
#$N25[1] 13023
#$N30[1] 12324
#$N35[1] 11682
#$N40[1] 11035
#$N45[1] 10363
#$N50[1] 9681
#$N55[1] 8992
#$N60[1] 8312
#$N65[1] 7638
#$N70[1] 6958
#$N75[1] 6269
#$N80[1] 5551
#$N85[1] 4756
#$N90[1] 3808
#$N95[1] 2613
