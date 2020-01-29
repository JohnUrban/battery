#!/bin/bash

echo BAM $BAM
echo BEDPE $BEDPE
echo MINSUPPORT ${MINSUPPORT}
echo OUTPRE $OUTPRE
echo

sniffles -m $BAM -b $BEDPE.bedpe --min_support $MINSUPPORT


## SVs from Sniffles
## NUM
grep -c -v ^# $BEDPE.bedpe > numsv

## SUM PREDICTED LENGTHS
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '{s+=$NF}END{print s}' > sumsv


## TYPES
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="DEL" {s+=1}END{print s}' > numdel
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="DUP" {s+=1}END{print s}' > numdup
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="INS" {s+=1}END{print s}' > numins
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="INV" || $11=="DEL/INV" || $11=="INVDUP" || $11=="INV/INVDUP" {s+=1}END{print s}' > numinv
grep -v ^# $BEDPE.bedpe | awk -v "s=0" '$11=="TRA" {s+=1}END{print s}' > numtra


## SV TYPES I ENCOUNTERED
#  DEL
#  DEL/INV
#  DUP
#  INS
#  INV
#  INVDUP
#  INV/INVDUP
#  TRA
