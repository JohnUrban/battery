#!/bin/bash

G=../../asms/refs/01-ecoligenome/01-ecoli_k12_mg1655_NC_000913.3.fasta

faSize -detailed $G > ecoli.genome


randomBed -l 1000 -n 25 -g ecoli.genome | fastaFromBed -bed - -fi $G > knownseqs.fa

