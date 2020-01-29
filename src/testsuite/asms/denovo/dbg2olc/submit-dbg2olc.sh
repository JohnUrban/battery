#!/bin/bash
###########################


#bash fragmentgenome.sh 


PB=`readlink -f ../../../data/sim_pb/pacbio.s.fasta`
PBONT=`readlink -f ../../../data/sim_pbont/pbont.fasta`
PLATANUS=`readlink -f inputasm.fasta`


sbatch -J dbg2olc_pb -c 16 --mem=100g --time=24:00:00 -o dbg2olc-pb-slurm-%A.out --export=PB=${PB},PLATANUS=${PLATANUS},DIR=pb dbg2olc.sh

sbatch -J dbg2olc_pbont -c 16 --mem=100g --time=24:00:00 -o dbg2olc-pbont-slurm-%A.out --export=PB=${PBONT},PLATANUS=${PLATANUS},DIR=pbont dbg2olc.sh



