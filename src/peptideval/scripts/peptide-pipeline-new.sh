#!/bin/bash
###########################

## $PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF $QUERYDIR $PRE $NJOBS

if [ $# -eq 0 ]; then echo "
Arg1 = scripts dir
Arg2 = config file
Arg3 = Logical(true/false) should directories be cleaned up...
Arg4 = QOS
Arg5 = /Path/To/Reference.fasta
Arg6 = Query Dir
Arg7 = Prefix to fasta files in query dir
Arg8 = Num Jobs -- equal to num files in query dir
Arg9 = Job Prefix
"; exit; fi

MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5
QUERYDIR=$6
PRE=$7
NJOBS=$8
JOBPRE=$9

BASE=`basename $ASM .fasta`

source $CONFIG

#Dir exist?
if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

if [ ! -d $BLASTOUTDIR ]; then
  mkdir $BLASTOUTDIR
fi
BLASTDIR=`readlink -f $MAIN`/$BLASTOUTDIR

if $CLEAN; then touch clean.txt; fi

##############################################################################
## MAKE BLAST DB FROM ASM
##############################################################################
D=blastdb
if $MAKEBLASTDB; then
 if [ ! -d $D ]; then mkdir $D; fi
 cd $D
 MAKEDONE=`sbatch -J ${JOBPRE}${BASE}_makeblastdb -o ${OUT}/makeblastdb.slurm.%A.out --mem=$MMEM --time=$MTIME -c $MTHREADS --account=${QOS} --export=ALL,ASM=${ASM},DBTYPE=nucl,OUT=asm $SCRIPTS/makeblastdb.sh | awk '{print $4}'`
 cd ../
fi
BDB=$(echo `readlink -f ${MAIN}/${D}`/asm)


##############################################################################
## BLAST
##############################################################################
TASK=blastn
BLASTDONE=`sbatch --dependency=afterok:${MAKEDONE} -a 1-$NJOBS -J ${JOBPRE}${BASE}_tblastn -o ${OUT}/tblastn.slurm.%A_%a.out --mem=$BMEM --time=$BTIME -c $BTHREADS --account=${QOS} \
   --export=ALL,QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},P=${BTHREADS},BDB=${BDB},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ} \
   ${SCRIPTS}/tblastn-array.sh | awk '{print $4}'`



##############################################################################
## FOLLOW UP 1
##############################################################################
FOLLOWUPNUM=1
FOLLOWDONE=`sbatch --dependency=afterany:${BLASTDONE} -J ${JOBPRE}${BASE}_tblastn_followup_${FOLLOWUPNUM} \
   -o ${OUT}/follow_up_${FOLLOWUPNUM}.slurm.%A.out \
   --mem=2g --time=6:00:00 -c 2 --account=${QOS} \
   --export=ALL,BASE=${BASE},NJOBS=${NJOBS},SLURMOUTDIR=${OUT},SLURMPRE=tblastn.slurm,FOLLOWUPNUM=${FOLLOWUPNUM},BMEM=${BMEM},BTIME=${BTIME},BTHREADS=${BTHREADS},QOS=${QOS},SCRIPTS=${SCRIPTS},QUERYDIR=${QUERYDIR},PRE=${PRE},BLASTDIR=${BLASTDIR},BDB=${BDB},EVAL=${EVAL},WORDSIZE=${WORDSIZE},CULL=${CULL},MAXTARGSEQ=${MAXTARGSEQ},JOBPRE=${JOBPRE} \
   ${SCRIPTS}/followup.sh | awk '{print $4}'`

