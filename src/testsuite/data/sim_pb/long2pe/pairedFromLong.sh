#!/bin/bash

FQ=../pacbio.s.fastq
NEW=pacbio
PRE=${NEW}-tmp

longRead2PairedReads.py -f $FQ --fq --ofq -l 1500 -r 500 -s 2000 -o $PRE
longRead2PairedReads.py -f $FQ --fq --ofq -l 1500 -r 500 -s 2000 -o $PRE.offset -S 500


cat ${PRE}-1.fastq ${PRE}.offset-1.fastq > ${NEW}-1.fastq
cat ${PRE}-2.fastq ${PRE}.offset-2.fastq > ${NEW}-2.fastq
cat ${PRE}.stats.txt ${PRE}.offset.stats.txt > ${NEW}.stats.txt
rm ${PRE}*


