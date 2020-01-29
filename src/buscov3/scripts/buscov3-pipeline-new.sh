#!/bin/bash
###########################
## Only differs from buscov3-pipeline.sh in changing --qos=$QOS to --account=$QOS
##$PIPELINE $SCRIPTS $CONFIG $CLEAN $QOS $REF

if [ $# -eq 0 ]; then echo "
Arg1 = SCRIPTS
Arg2 = CONFIG
Arg3 = CLEAN
Arg4 = QOS
Arg5 = REF
"; exit; fi


MAIN=$PWD

SCRIPTS=$1
CONFIG=$2
CLEAN=$3
QOS=$4
ASM=$5
BASE=`basename $ASM .fasta`

source $CONFIG

if [ ! -d $SLURMOUTDIR ]; then mkdir $SLURMOUTDIR; fi
OUT=`readlink -f $MAIN`/$SLURMOUTDIR

### PIPELINE

##############################################################################
## BUSCO V3 LOOP
##############################################################################

for DIR in ${LINEAGE_OPTIONS}; do
  LINEAGE=${LINEAGEBASE}/${DIR}_odb9

  #for debugging: for VAR in L LINEAGE RUN; do echo $VAR ${!VAR}; done

  if [ ! -d $DIR ]; then mkdir $DIR; fi
  cd $DIR
  DONE=`sbatch -J ${BASE}_buscov3_${DIR} -o ${OUT}/${DIR}.slurm.%A.out --mem=$BV3_MEM --time=$BV3_TIME -c $BV3_THREADS --account=$QOS --export=ALL,FASTA=${ASM},OUT=${BASE},CPU=${BV3_THREADS},LINEAGE=${LINEAGE},BV3_MODE=${BV3_MODE},REGIONLIMIT=${REGIONLIMIT} ${SCRIPTS}/buscov3.eval.sh | awk '{print $4}'`
  if $CLEAN; then
    MOP=`sbatch --dependency=afterok:${DONE} -J ${BASE}_buscov3_${DIR}_clean -o ${OUT}/clean_${DIR}.slurm.%A.out --mem=1g --time=1:00:00 -c 1 --account=$QOS ${SCRIPTS}/buscov3.clean.sh | awk '{print $4}'` 
  cd ../
 fi
done


