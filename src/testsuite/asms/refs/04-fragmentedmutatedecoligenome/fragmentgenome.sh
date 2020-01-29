#!/bin/bash

python2.7 -c "
from Bio import SeqIO
f = '../03-mutatedecoligenome/03-ecoli.mutated.fasta'
fa = SeqIO.read(f,'fasta')

orig = fa.id
fa.description = ''
step=500000
with open('04-ecoli-fragmented-mutated-genome.fasta','w') as out:
  for i in range(0,len(fa)-10000,step):
    fa.id = orig + '_' + str(i)
    SeqIO.write(fa[i:i+step], out, 'fasta')
  
"
