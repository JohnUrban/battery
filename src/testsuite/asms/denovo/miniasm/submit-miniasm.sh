#!/bin/bash
###########################

PB=`readlink -f ../../../data/sim_pb/pacbio.s.fastq`
PBONT=`readlink -f ../../../data/sim_pbont/pbont.fastq`

T=8

DIR=pb
sbatch -J miniasm_pb -c 8 --mem=64g --time=12:00:00 -o miniasm-pb-slurm-%A.out --export=READS=${PB},DIR=${DIR},T=${T} miniasm.sh


DIR=pbont
sbatch -J miniasm_pbont -c 8 --mem=64g --time=12:00:00 -o miniasm-pbont-slurm-%A.out --export=READS=${PBONT},DIR=${DIR},T=${T} miniasm.sh


