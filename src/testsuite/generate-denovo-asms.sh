#!/bin/bash
set -e

if [ $# -eq 0 ]; then echo "Need first/only argument: /path/to/battery_paths.cfg"; exit; fi
source ${1}

## clear old files if present
bash reset.sh

bash third_party/check-python.sh
bash third_party/check-various.sh

##################################
### DE NOVO ASMS #################
##################################

## go into denovo
SBATCH=`command -v sbatch`
if [ ! -z $SBATCH ]; then
  echo "...Starting De Novo Assemblies"
  #enter denovo
  cd asms/denovo
  #process abruijn
  cd abruijn
  echo "...Starting De Novo Assemblies: abruijn"
  bash submit-abruijn.sh
  cd ../
  # process canu
  cd canu
  echo "...Starting De Novo Assemblies: canu"
  bash submit-canu.sh 
  cd ../
  #process dbg2olc
  cd dbg2olc
  echo "...Starting De Novo Assemblies: dbg2olc"
  bash fragmentgenome.sh
  bash submit-dbg2olc.sh
  cd ../
  # process falcon
  cd falcon
  echo "...Starting De Novo Assemblies: falcon"
  bash submit-falcon.sh
  cd ../
  #process miniasm
  cd miniasm
  echo "...Starting De Novo Assemblies: miniasm"
  bash submit-miniasm.sh
  cd ../  
  #leave asms/denovo - back to base
  cd ../../
fi

