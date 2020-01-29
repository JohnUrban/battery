#!/bin/bash


COROUTCOV=500

module load java/8u66

T=24:00:00
QOS=epscor-condo

G=4.6m

/gpfs_home/jurban/software/canu/canu/Linux-amd64/bin/canu \
 -p $NAME -d $NAME \
 genomeSize=$G \
 -pacbio-raw $PB \
 "gridOptions=--time $T --qos=$QOS" \
 oeaMemory=8 \
 corOutCoverage=$COROUTCOV \
 corMemory=30
