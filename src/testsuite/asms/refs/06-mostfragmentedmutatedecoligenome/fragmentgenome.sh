#!/bin/bash

python2.7 -c "
from Bio import SeqIO
#f = '../mutatedecoligenome/03-ecoli.mutated.fasta'
#f = '../fragmentedmutatedecoligenome/04-ecoli-fragmented-mutated-genome.fasta'
f = '../05-mutatefragmentedgenomemore/05-ecoli.mutate-frag-mutate.fasta'

step=250000
with open('06-ecoli-most-fragmented-and-mutated-genome.fasta','w') as out:
  for fa in SeqIO.parse(f, 'fasta'):
    orig = fa.id
    fa.description = ''
    for i in range(0,len(fa)-10000,step):
      fa.id = orig + '_' + str(i)
      SeqIO.write(fa[i:i+step], out, 'fasta')
  
"
