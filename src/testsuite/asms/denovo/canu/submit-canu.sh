#!/bin/bash
###########################


PB=`readlink -f ../../../data/sim_pb/pacbio.s.fasta`
PBONT=`readlink -f ../../../data/sim_pbont/pbont.fasta`
COROUTCOV=500
module load java/8u66
T=24:00:00
QOS=epscor-condo
G=4.6m




NAME=pb
/gpfs_home/jurban/software/canu/canu/Linux-amd64/bin/canu \
 -p $NAME -d $NAME \
 genomeSize=$G \
 -nanopore-raw $PB \
 "gridOptions=--time $T --qos=$QOS" \
 corOutCoverage=$COROUTCOV 


NAME=pbont
/gpfs_home/jurban/software/canu/canu/Linux-amd64/bin/canu \
 -p $NAME -d $NAME \
 genomeSize=$G \
 -nanopore-raw $PBONT \
 "gridOptions=--time $T --qos=$QOS" \
 corOutCoverage=$COROUTCOV 


# oeaMemory=8 \
# corMemory=30
