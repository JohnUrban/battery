#!/bin/bash
###########################

PB=`readlink -f ../../../data/sim_pb/pacbio.s.fasta`
PBONT=`readlink -f ../../../data/sim_pbont/pbont.fasta`

T=32


COV=50
DIR=pb
sbatch -J abruijn_pb -c 32 --mem=200g --time=48:00:00 -o abruijn-pb-slurm-%A.out --export=READS=${PB},DIR=${DIR},T=${T},COV=${COV} abruijn.sh


COV=60
DIR=pbont
sbatch -J abruijn_pbont -c 32 --mem=200g --time=48:00:00 -o abruijn-pbont-slurm-%A.out --export=READS=${PBONT},DIR=${DIR},T=${T},COV=${COV} abruijn.sh


