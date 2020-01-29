#!/bin/bash
set -e


if [ $# -eq 0 ]; then echo "Need first/only argument: /path/to/battery_paths.cfg"; exit; fi
source ${1}

## DOES NOT DO THE DENOVO ASM LAUNCHES

## clear old files if present
bash reset.sh

bash third_party/check-python.sh
bash third_party/check-various.sh


##################################
### ASSEMBLIES ###################
##################################

echo "...Generating Assemblies"

## go into asm folder - references
cd asms/refs/

## generate needed files

echo "...Generating Assemblies: snpd"
cd 02-snpdecoligenome/
bash snpthegenome.sh
cd ../

echo "...Generating Assemblies: snpd and mutated"
cd 03-mutatedecoligenome/
bash modifygenome.sh
cd ../

echo "...Generating Assemblies: snpd, mutated, and fragmented"
cd 04-fragmentedmutatedecoligenome
bash fragmentgenome.sh
cd ../

echo "...Generating Assemblies: snpd, mutated, fragmented, and mutated more"
cd 05-mutatefragmentedgenomemore
bash modifygenome.sh
cd ../

echo "...Generating Assemblies: most fragmented and mutated"
cd 06-mostfragmentedmutatedecoligenome/
bash fragmentgenome.sh
cd ../

## leave asm folder - references
cd ../../



##################################
### DATA #########################
##################################

echo "...Generating Data"

## enter data folder
cd data/

cd sim_ilmn_pe
echo "...Generating Data: illumina"
bash makereads.sh
cd ../

cd sim_ont
echo "...Generating Data: ont"
bash makereads.sh
cd long2pe
echo "...Generating Data: ont long2pe"
bash pairedFromLong.sh
cd ../../

cd sim_pb
echo "...Generating Data: pb"
bash makereads.sh
cd long2pe
echo "...Generating Data: pb long2pe"
bash pairedFromLong.sh
cd ../../

cd sim_pbont
echo "...Generating Data: pbont"
bash get.sh
cd ../

cd ecoliknownseqs
echo "...Generating Data: knownseqs"
bash getknowseqs.sh 
cd ../

cd sim_rna_seq
echo "...Generating Data: rnaseq"
bash makereads.sh
cd ../


## leave data folder
cd ../




##################################
### DATA - BNG ###################
##################################
## Takes a long time -- so I let the denovo asms start before hand

cd data/sim_bng
echo "...Generating Data: bng"
bash makereads.sh
cd ../../



	
##################################
### EVALS ########################
##################################
bash generate-evals-dir.sh 

