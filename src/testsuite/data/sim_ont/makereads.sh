#!/bin/bash

G=../../asms/refs/01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta
O=ont
L=8000
C=10
#fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

LENS="59190 39784 31338 25718 21477 18036 15198 12964 11261 9934 8869 7958 7143 6360 5591 4782 3928 3017 2027 1000"
for L in $LENS; do
  EM=0.05
  EI=0.05
  ED=0.05
  C=0.02
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.06
  EI=0.07
  ED=0.07
  C=0.03
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.07
  EI=0.09
  ED=0.09
  C=0.05
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.1
  EI=0.1
  ED=0.1
  C=0.11
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.11
  EI=0.12
  ED=0.12
  C=0.11
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.13
  EI=0.14
  ED=0.13
  C=0.1
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se

  EM=0.15
  EI=0.15
  ED=0.15
  C=0.1
  O=tmp-${L}-${C}-${EM}-${EI}-${ED}
  fastqSimulate -f $G -o $O -l $L -x $C -em $EM -ei $EI -ed $ED -se
done 2> errfile
cat tmp*fastq > ont.s.fastq
rm tmp*fastq

#$N5[1] 59190
#$N10[1] 39784
#$N15[1] 31338
#$N20[1] 25718
#$N25[1] 21477
#$N30[1] 18036
#$N35[1] 15198
#$N40[1] 12964
#$N45[1] 11261
#$N50[1] 9934
#$N55[1] 8869
#$N60[1] 7958
#$N65[1] 7143
#$N70[1] 6360
#$N75[1] 5591
#$N80[1] 4782
#$N85[1] 3928
#$N90[1] 3017
#$N95[1] 2027
