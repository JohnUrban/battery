#!/bin/bash

# specify paths to lap read sample (LR1,LR2) and all reads (R1,R2)-- give dummy answers if will not be using (that will serve as place-holders)
LR1=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.1.fastq
LR2=/users/jurban/data/scratch/lap/sample-1.5m/downsampled.2.fastq
R1=~/data/scratch/male-ilmn/data/ilmnraw/R1.fastq
R2=~/data/scratch/male-ilmn/data/ilmnraw/R2.fastq


## What programs to use? FILL IN BELOW
ALL=eval.cfg
OnlyAle=eval.aleonly.cfg
OnlyBusco=eval.buscoOnly.cfg
OnlyLap=eval.laponly.cfg
OnlyReapr=eval.reapronly.cfg
OnlyReaprNoClean=eval.reapronly.noclean.cfg
OnlyReaprNoCleanAggressive=eval.reapronly.noclean.aggressive.cfg

## FILL IN WITH CORRECT VARIABLE
EvalThese=$ALL



## May need to adjust the following
ASMDIR=asms
SCRIPTS=/gpfs_home/jurban/software/sciaratools/sciara-project-tools/slurmgear/shortreadeval/scripts/
EVAL=${SCRIPTS}/eval.ARGS.sh
CONFIG=${SCRIPTS}/configs/${EvalThese}


# eval.ARGS.sh
#Arg1=/Path/To/Reference.fasta
#Arg2=QOS
#Arg3=scripts dir path
#Arg4=lap reads 1
#Arg5=lap reads 2
#Arg6=all reads 1
#Arg7=all reads 2
#Arg8=ConfigFile


i=0

for f in ${ASMDIR}/*fasta; do 
  i=$(( $i+1 ))
  if [ $i -eq 9 ]; then QOS=biomed-sb-condo; i=0; else QOS=epscor-condo; fi
  ref=`readlink -f $f` 
  b=`basename $f .fasta`; 
  echo $b; 
  mkdir $b; 
  cd $b;
  $EVAL $ref $QOS $SCRIPTS $LR1 $LR2 $R1 $R2 $CONFIG
  cd ../
done 

