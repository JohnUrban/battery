#!/bin/bash

## note that error is not simulated here
## RMAP format name, bp len, nfrags, fraglen1.....fraglenN

G=../../asms/refs/01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta
N=4000
L=200000
##REC_SEQ=GGATCC # recognition sequence for BamHI
REC_SEQ=CACGAG # recognition sequence for BssSI (Sciara data)

MIN_FRAG_SIZE=1000 # bp units


# MAKE GENOME FILE
faSize -detailed $G > ecoli.genome

## GET FASTA SEQS TO BE CONVERTED TO RMAPS
bedtools random -g ecoli.genome -l $L -n $N | fastaFromBed -fi $G -bed - | fasta_name_changer.py -f - --replace rmap -n  > ecoli.rmaps.fasta


# CONVERT FASTA SEQS TO RMAPS
make_insilico_map -o ecoli ecoli.rmaps.fasta $REC_SEQ 2> conversion.output.txt
smooth_maps_file -m $MIN_FRAG_SIZE ecoli.maps > ecoli.smoothed.maps


echo  ${PWD}/ecoli.smoothed.maps > maps.fofn
