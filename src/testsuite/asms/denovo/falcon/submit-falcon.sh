#!/bin/bash
###########################


PB=`readlink -f ../../../data/sim_pb/pacbio.s.fasta`
PBONT=`readlink -f ../../../data/sim_pbont/pbont.fasta`



sbatch -o slurm-pb-%A.out -J falcon_pb -c 4 --mem=12g --time=48:00:00 --export=READS=${PB},DIR=pb falcon.sh
sbatch -o slurm-pbont-%A.out  -J falcon_pb -c 4 --mem=12g --time=48:00:00 --export=READS=${PBONT},DIR=pbont falcon.sh
