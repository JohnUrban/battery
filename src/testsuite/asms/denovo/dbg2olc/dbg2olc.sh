#!/bin/bash

## import PLATANUS, PB, DIR

mkdir $DIR
cd $DIR

echo dbg2olc; mkdir dbg2olc; cd dbg2olc
ln -s /users/jurban/software/dbg2olc/DBG2OLC/compiled/* .
./DBG2OLC k 17 KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.002 LD1 0 MinLen 200 Contigs $PLATANUS RemoveChimera 1 f $PB
cd ../; echo; echo; echo; echo

echo cat
cat $PLATANUS $PB > ctg_pb.fasta
echo; echo; echo; echo

echo pbdagcon_consensus; mkdir pbdagcon_consensus; cd pbdagcon_consensus
module load blasr/2015Oct22-8cc8621
PATH=~/software/pbdagcon/pbdagcon/src/cpp/:~/software/pbdagcon/pbdagcon/src/:~/software/dbg2olc/johnfork/DBG2OLC/utility/:$PATH
ln -s /users/jurban/software/dbg2olc/dbg2olcscripts-fromJJemerson-Mahul-etal/* .
mkdir consensus_dir
split_and_run_pbdagcon.path.sh ../dbg2olc/backbone_raw.fasta ../dbg2olc/DBG2OLC_Consensus_info.txt ../ctg_pb.fasta consensus_dir > consensus_log.txt
